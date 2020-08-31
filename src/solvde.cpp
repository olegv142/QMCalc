#include "errs.h"
#include "nrutil.h"
#include "nrtools.h"

#include <stdio.h>
#include <math.h>

/* Driver routine for solution of two point boundary value problems by relaxation. ITMAX is the
maximum number of iterations. CONV is the convergence criterion. SLOWC controls
the fraction of corrections actually used after each iteration. SCALV[1..ne] contains typical
sizes for each dependent variable, used to weight errors. INDEXV[1..ne] lists the column
ordering of variables used to construct the matrix S[1..ne][1..2*ne+1] of derivatives. ( The
NB boundary conditions at the first mesh point must contain some dependence on the first NB
variables listed in INDEXV ). The problem involves NE equations for NE adjustable dependent
variables at each point. At the first mesh point there are NB boundary conditions. There are a
total of M mesh points. Y[1..ne][1..m] is the two-dimensional array that contains the initial
guess for all the dependent variables at each mesh point. On each iteration, it is updated by the
calculated correction. The arrays C[1..ne][1..ne-nb+1][1..m+1] and s supply dummy
storage used by the relaxation code.
  The only information passed from difeq to solvde is the matrix of derivatives
S[1..ne][1..2*ne+1]; all other arguments are input to difeq and should not be altered.
K indicates the current mesh point, or block number. If K = 1 or K > M, the block
involves the boundary conditions at the first or final points; otherwise the block
acts on FDEs coupling variables at point (K-1,K).
  In matrix S[i][j] rows I label equations, columns J refer to
derivatives with respect to dependent variables in the solution. Each equation will
depend on the NE dependent variables at either one or two points. Thus j runs from 1
to 2 * NE. The column ordering for dependent variables at each point must agree
with the list supplied in INDEXV[j]. Thus, for a block not at a boundary, the first
column multiplies dY(l=indexv[1],k-1), and the column NE+1 multiplies
dY(l=indexv[1],k). The difference equations must be stored in column [2*ne+1].
Boundary conditions must always be specifyed in S[i][ne..2*ne] part of the S matrix.
Left boundary conditions must be specified in last NB rows while right
boundary - in first NE-NB rows. The rest part of S in such cases may be arbitrary.
The NB left boundary conditions must contain some dependence on the first NB 
dependent variables. If not, then the total matrix will appear to be singular
and the method will fail. */

DefineErrorMess( SOLVDE_IT, "Too many iterations in solvde()" );

int solvde( DiffEqCb difeq, int itmax, float conv, float slowc, float scalv[], int indexv[],
	int ne, int nb, int m, float **y, float ***c, float **s, FILE *debug, void* ctx )
{
	void bksub( int ne, int nb, int jf, int k1, int k2, float ***c );
	void pinvs( int ie1, int ie2, int je1, int jsf, int jc1, int k, float ***c, float **s );
	void red( int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
			  int ic1, int jc1, int jcf, int kc, float ***c, float **s );

	int it, j, jv, k, km, kp;
	float err, errj, fac, vmax, vz;

	int   *kmax = ivector( 1, ne );
	float *ermax = vector( 1, ne );

	int nvars = ne * m;
	int nb1 = nb + 1;
	int ne1 = ne + 1;
	int neb = ne + nb;
	int neb1 = neb + 1;
	int n2e = ne + ne;
	int n2e1 = n2e + 1;
	int nc = ne - nb;
	int nc1 = nc + 1;

	for( it = 1 ; it <= itmax ; it++ ) {    /* Primary iteration loop */
		try {
	                k = 1;                             /* Boundary conditions at first point */
			difeq( k, indexv, s, y, ctx );
			pinvs( nc1, ne, ne1, n2e1, 1, 1, c, s );
	                for( k = 2 ; k <= m ; k++ ) { /* Finite difference equations at all point pairs */
		                kp = k - 1;
				difeq( k, indexv, s, y, ctx );
				red( 1, ne, 1, nb, nb1, ne, n2e1, nc1, 1, nc1, kp, c, s );
				pinvs( 1, ne, nb1, n2e1, 1, k, c, s );
	                }
		        k = m + 1;  /* Final boundary conditions */
			difeq( k, indexv, s, y, ctx );
			red( 1, nc, ne1, neb, neb1, n2e, n2e1, nc1, 1, nc1, m, c, s );
			pinvs( 1, nc, neb1, n2e1, nc1, k, c, s );
		} _catch {
			if( lastError() != INTERNAL_ERROR ) {
				free_vector( ermax, 1, ne );
				free_ivector( kmax, 1, ne );
			}
			throw;
		}
		bksub( ne, nb, nc1, 1, m, c );  /* Backsubstitution */
		err = 0;
		for( j = 1 ; j <= ne ; j++ ) { /* Convergence check, accumulate error */
			jv = indexv[j];
			errj = vmax = 0;
			km = 0;
			for( k = 1 ; k <= m ; k++ ) { /* Find point with lagest error, */
				vz = fabs( c[jv][1][k] );   /* for each dependent variable */
				if( vz > vmax ) {
					vmax = vz;
					km = k;
				}
				errj += vz;
			}
			err += errj / scalv[j];   /* Note weighting for each dependent variable */
			ermax[j] = c[jv][1][km] / scalv[j];
			kmax[j] = km;
		}
		err /= nvars;
		fac = ( err > slowc ? slowc / err : 1 );
		for( j = 1 ; j <= ne ; j++ ) { /* Apply corrections */
			jv = indexv[j];
			for( k = 1 ; k <= m ; k++ )
			y[j][k] -= fac * c[jv][1][k];
		}
		if( debug ) {
			fprintf( debug, "*%8s %9s %9s\n", "Iter.", "Error", "FAC" ); /* Summary */
			fprintf( debug, "*%6d %12.6f %11.6f\n", it, err, fac );
		}
		if( err < conv ) {
			free_vector( ermax, 1, ne );
			free_ivector( kmax, 1, ne );
			return it;
		}
	}
	free_vector( ermax, 1, ne );
	free_ivector( kmax, 1, ne );
	Signal( SOLVDE_IT );  /* Convergence failed */
	return it;
}

void bksub( int ne, int nb, int jf, int k1, int k2, float ***c )
/* Backsubstitution, used internally by solvde */
{
	int nbf, im, kp, k, j, i;
	float xx;

	nbf = ne - nb;
	im = 1;
	for( k = k2 ; k >= k1 ; k-- ) {  /* Use recurrence relations to eliminate */
		if( k == k1 ) im = nbf + 1;  /* remaining dependences */
		kp = k + 1;
		for( j = 1 ; j <= nbf ; j++ ) {
			xx = c[j][jf][kp];
			for( i = im ; i <= ne ; i++ )
				c[i][jf][k] -= c[i][j][k] * xx;
		}
	}
	for( k = k1 ; k <= k2 ; k++ ) { /* Reorder corrections to be in column 1 */
		kp = k + 1;
		for( i = 1 ; i <= nb ; i++ ) c[i][1][k] = c[i+nbf][jf][k];
		for( i = 1 ; i <= nbf ; i++ ) c[i+nb][1][k] = c[i][jf][kp];
	}
}

void pinvs( int ie1, int ie2, int je1, int jsf, int jc1, int k, float ***c, float **s )
/* Diagonalize the square subsection of the s matrix, and store the recursion coefficients in C;
used internally by solvde. */
{
	int js1, jpiv, jp, je2, jcoff, j, irow, ipiv, id, icoff, i, *indxr;
	float pivinv, piv, big, *pscl;

	indxr = ivector( ie1, ie2 );
	pscl = vector( ie1, ie2 );
	je2 = je1 + ie2 - ie1;
	js1 = je2 + 1;
	for( i = ie1 ; i <= ie2 ; i++ ) {  /* Implicit pivoting */
		big = 0;
		float tmp;
		float *ptr = s[i] + je1;
		float *end = s[i] + je2;
		while( ptr <= end )
			if( ( tmp = fabs( *ptr++ ) ) > big ) big = tmp;
		if( big == 0 ) {
			free_vector( pscl, ie1, ie2 );
			free_ivector( indxr, ie1, ie2 );
			nrerror( "Singular matrix - row all 0 in pinvs()" );
		}
		pscl[i] = 1 / big;
		indxr[i] = 0;
	}
	for( id = ie1 ; id <= ie2 ; id++ ) {
		piv = 0;
		float **rptr = s + ie1;
		for( i = ie1 ; i <= ie2 ; i++ ) {  /* Find pivot element */
			if( indxr[i] == 0 ) {
				big = 0;
				float *cptr = *rptr + je1;
				for( j = je1 ; j <= je2 ; j++ ) {
					float tmp = fabs( *cptr );
					if( tmp > big ) {
						jp = j;
						big = tmp;
					}
					cptr++;
				}
				float tmp = big * pscl[i];
				if( tmp > piv ) {
					ipiv = i;
					jpiv = jp;
					piv = tmp;
				}
			}
			rptr++;
		}
		if( s[ipiv][jpiv] == 0 ) {
			free_vector( pscl, ie1, ie2 );
			free_ivector( indxr, ie1, ie2 );
			nrerror( "Singular matrix in pinvs()" );
		}
		indxr[ipiv] = jpiv; /* In place reduction. Save column ordering. */
		pivinv = 1 / s[ipiv][jpiv];
		float *ptr = s[ipiv] + je1;
		float *end = s[ipiv] + jsf;
		while( ptr <= end ) *ptr++ *= pivinv; /* Normalize pivot row */
		s[ipiv][jpiv] = 1;
		for( i = ie1 ; i <= ie2 ; i++ ) {  /* Reduce nonpivot elements in column */
			if( indxr[i] != jpiv ) {
				float *si = s[i];
				float &dum = si[jpiv];
				if( dum ) {
					float tmp = dum;
					float *aptr = si + je1;
					float *bptr = s[ipiv] + je1;
					float *cptr = si + jsf;
					while( aptr <= cptr )
						*aptr++ -= tmp * *bptr++;
					dum = 0;
				}
			}
		}

	}
	jcoff = jc1 - js1;  /* Sort and store unreduced coefficients */
	icoff = ie1 - je1;
	for( i = ie1 ; i <= ie2 ; i++ ) {
		irow = indxr[i] + icoff;
		float  *sptr = s[i] + js1;
		float **cptr = c[irow] + js1 + jcoff;
		float **eptr = c[irow] + jsf + jcoff;
		while( cptr <= eptr ) (*cptr++)[k] = *sptr++;
	}
	free_vector( pscl, ie1, ie2 );
	free_ivector( indxr, ie1, ie2 );
}

void red( int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
			  int ic1, int jc1, int jcf, int kc, float ***c, float **s )
/* Reduce columns jz1..jz2 of the S matrix, using previous results as stored in C matrix. Only
columns jm1..jm2, jmf are affected by the prior results. red is used internally by solvde. */
{
	int loff, l, j, ic, i;
	float vx;

	loff = jc1 - jm1;
	ic = ic1;
	for( j = jz1 ; j <= jz2 ; j++ ) { /* Loop over columns to be zeroed */
		for( l = jm1 ; l <= jm2 ; l++ ) { /* Loop over columns altered */
			vx = c[ic][l+loff][kc];
			for( i = iz1 ; i <= iz2 ; i++ ) s[i][l] -= s[i][j] * vx;  /* Loop over rows */
		}
		vx = c[ic][jcf][kc];
		for( i = iz1 ; i <= iz2 ; i++ ) s[i][jmf] -= s[i][j] * vx;  /* Plus final element */
		ic += 1;
	}
}
