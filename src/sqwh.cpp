#include "sqwh.h"
#include "nrutil.h"
#include "nrtools.h"
#include "express.h"
#include "misc.h"
#include <fstream>  

#define _USE_MATH_DEFINES
#include <math.h>

#define ITMAX 1000
#define CONV0 .1f
#define SLOWC CONV0

// The number of equations
#define NE 10
/* The number of equations matches the number of independent variables:
  [1]     The wavefunction density integral (0 at left bound, 1 at right bound)
  [2..5]  Wavefunction spin components
  [6..9]  The wavefunction 'derivative' spin components (actually the result of applying speed operator to wavefunction, see description in doc)
  [10]    Energy eigenvalue
*/

// The number of boundary conditions at left boundary
#define NB 5

DefineErrorMess( SQWH_OUT, "Can't write output file" );
DefineErrorMess( SQWH_QNCH, "Wavefunction quench" );

void GetSQWHParams(struct SQWHParams& p, TExpression const& e)
{
	p.M        = (unsigned)e.Get("grid");
	p.NL       = (unsigned)e.Get("NL");
	p.subband  = (unsigned)e.Get("subband");
	p.psteps   = (unsigned)e.Get("psteps");
	p.W        = (float   )e.Get("W");
	p.hb0      = (float   )e.Get("hb0");
	p.hb1      = (float   )e.Get("hb1");
	p.Bstep    = (float   )e.Get("Bstep");
	p.Bsteps   = (unsigned)e.Get("Bsteps");
	p.Eo       = (float   )e.Get("Eo");
	p.Bo       = (float   )e.Get("Bo");
	p.g1       = (float   )e.Get("g1");
	p.g2       = (float   )e.Get("g2");
	p.g3       = (float   )e.Get("g3");
	p.K        = (float   )e.Get("K");
	p.prec     = (float   )e.Get("prec");
}

SQWHSolver::SQWHSolver(struct SQWHParams const& params)
	: p(params)
{
	y = matrix( 1, NE, 1, p.M );
	s = matrix( 1, NE, 1, 2 * NE + 1 );
	c = f3tensor( 1, NE, 1, NE - NB + 1, 1, p.M + 1 );

	pot    = vector( 1, p.M );
	scalv  = vector( 1, NE );
	indexv = ivector( 1, NE );

	mh = 1 / (p.g1 / 2 - p.g2);
	ml = 1 / (p.g1 / 2 + p.g2);
	A = p.g1 + p.g2;
	B = p.g1 - p.g2;

	unsigned i;
	for( i = 1; i <= NE ; i++ ) {
		scalv[i] = 1;
		indexv[i] = i;
	}
	scalv[NE] = 1 / (ml + mh);

	float B2 = ( 1 - p.W ) / 2;
	for( i = 1 ; i <= p.M ;  i++ ) {
		float z = ( i - 1.f ) / ( p.M - 1.f );
		if( z < B2 ) pot[i] = p.hb0;
		else if( z > 1 - B2 ) pot[i] = p.hb1;
		else pot[i] = 0;
	}

	D = matrix( 1, NE, 1, NE );

	z[0] = matrix( 1, NE, 1, p.M );
	z[1] = matrix( 1, NE, 1, p.M );
	z[2] = matrix( 1, NE, 1, p.M );
	z[3] = matrix( 1, NE, 1, p.M );

	if (p.NL > 0) {
		e = f3tensor( 0, p.NL - 1, 0, p.Bsteps, 0, 3 );
		clear_f3tensor( e, 0, p.NL - 1, 0, p.Bsteps, 0, 3 );
	} else
		e = NULL;
}

SQWHSolver::~SQWHSolver()
{
	free_matrix( y, 1, NE, 1, p.M );
	free_matrix( s, 1, NE, 1, 2 * NE + 1 );
	free_f3tensor( c, 1, NE, 1, NE - NB + 1, 1, p.M + 1 );
	free_vector( pot, 1, p.M );
	free_vector( scalv, 1, NE );
	free_ivector( indexv, 1, NE );
	free_matrix( D, 1, NE, 1, NE );
	free_matrix( z[0], 1, NE, 1, p.M );
	free_matrix( z[1], 1, NE, 1, p.M );
	free_matrix( z[2], 1, NE, 1, p.M );
	free_matrix( z[3], 1, NE, 1, p.M );
	if (p.NL > 0)
		free_f3tensor( e, 0, p.NL - 1, 0, p.Bsteps, 0, 3 );
}

// Create initial guess for the particular spin s = 0..3
void SQWHSolver::init_guess(unsigned spin)
{
	float m = (spin == 0 || spin == 3) ? mh : ml;
	float pl = (float)M_PI * (p.subband + 1.f);
	float e = pl * pl / m;
	unsigned i;
	clear_matrix( y, 1, NE, 1, p.M );
	for( i = 1 ; i <= p.M ;  i++ ) {
		float z = ( i - 1.f ) / ( p.M - 1.f );
		float f = pl * z;
		y[1][i] = z - sin(2 * f) / (2 * pl);
		y[2 + spin][i] = sqrt(2.f) * sin(f);
		y[6 + spin][i] = sqrt(2.f) * pl * cos(f) / m;
		y[10][i] = e;
	}
	pscale = 0;
}

void SQWHSolver::init_level(unsigned lvl)
{
	n  = lvl;
	n0 = lvl + .5f;
	n1 = lvl > 0 ? lvl - 1 + .5f : 0;
	n2 = lvl > 1 ? lvl - 2 + .5f : 0;
	n3 = lvl > 2 ? lvl - 3 + .5f : 0;

	float R = sqrt(3.f) * (p.g2 + p.g3) / 2;
	float S = sqrt(3/2.f) * p.g3;
	M = n > 2 ? R * sqrt((n - 1.f)*(n - 2.f)) : 0;
	N = n > 1 ? R * sqrt(n*(n - 1.f)) : 0;
	P = n > 2 ? S * sqrt(n - 2.f) : 0;
	Q = n > 0 ? S * sqrt((float)n) : 0;
	P2 = P * P;
	Q2 = Q * Q;

	set_magnetic_field(0);
}

void SQWHSolver::set_magnetic_field(float field)
{
	H = field;
	sH = sqrt(field);
	init_derivative_matrix();
}

void SQWHSolver::Solve()
{
	solve_zero_field();
	for (unsigned l = 0; l < p.NL; ++l)
		solve_level(l);
}

void SQWHSolver::SaveResults(const char* basename) const
{
	std::string bname(basename);
	bname += (char)('0' + p.subband);

	save_wavefunction(z[0], bname + "_hm.dat");
	save_wavefunction(z[1], bname + "_lm.dat");
	save_wavefunction(z[2], bname + "_lp.dat");
	save_wavefunction(z[3], bname + "_hp.dat");

	bname += 'L';
	for (unsigned l = 0; l < p.NL; ++l) {
		save_level(e[l], bname + (char)('0' + l) + ".dat");
	}
}

// Count zeros for the particular spin component s = 0..3
unsigned SQWHSolver::count_zeros(unsigned spin, float precision) const
{
	return ::count_zeros(y[2+spin], 1, p.M, precision);
}

/*
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
and the method will fail.
*/
void SQWHSolver::eq(int k, int* idx, float **s, float **y) const
{
	float h = 1 / (p.M - 1.f), h2 = h / 2;
	clear_matrix( s, 1, NE, 1, 2 * NE + 1 );

#define EQ (2*NE+1)
	if( k == 1 ) {
		// The y[1..5] = 0 at left boundary
		s[6][EQ]  = y[1][1];
		s[7][EQ]  = y[2][1];
		s[8][EQ]  = y[3][1];
		s[9][EQ]  = y[4][1];
		s[10][EQ] = y[5][1];
		// Derivatives:
		s[6][NE+1] = s[7][NE+2] = s[8][NE+3] = s[9][NE+4] = s[10][NE+5] = 1;
	} else if( k > (int)p.M ) {
		// The y[1] = 1, y[2..5] = 0 at right boundary
		s[1][EQ]  = y[1][p.M] - 1;
		s[2][EQ]  = y[2][p.M];
		s[3][EQ]  = y[3][p.M];
		s[4][EQ]  = y[4][p.M];
		s[5][EQ]  = y[5][p.M];
		// Derivatives:
		s[1][NE+1] = s[2][NE+2] = s[3][NE+3] = s[4][NE+4] = s[5][NE+5] = 1;
	} else {
		// Internal pair of points
		float a[NE+1]; // 2 point average
		float d[NE+1]; // 2 point deltas
		int i, j;

		// Fill diagonal elements first
		for( i = 1 ; i <= NE ; i++ ) {
			a[i] = ( y[i][k] + y[i][k-1] ) / 2;
			d[i] = y[i][k] - y[i][k-1];
			s[i][i] = -1;
			s[i][NE+i] = 1;
		}

		// Fill equations
		s[1][EQ] = d[1] - h * (a[2]*a[2] + a[3]*a[3] + a[4]*a[4] + a[5]*a[5]);
		s[2][EQ] = d[2] - h * mh * (a[6] - sH * Q * a[3]);
		s[3][EQ] = d[3] - h * ml * (a[7] + sH * Q * a[2]);
		s[4][EQ] = d[4] - h * ml * (a[8] + sH * P * a[5]);
		s[5][EQ] = d[5] - h * mh * (a[9] - sH * P * a[4]);

		float phy = pscale * (pot[k] + pot[k-1]) / 2;
		float E0 = (A*n0 - ml*Q2 - 3*p.K/2)*H + phy - a[10];
		float E1 = (B*n1 - mh*Q2 -   p.K/2)*H + phy - a[10];
		float E2 = (B*n2 - mh*P2 +   p.K/2)*H + phy - a[10];
		float E3 = (A*n3 - ml*P2 + 3*p.K/2)*H + phy - a[10];
	
		s[6][EQ] = d[6] - h * (E0 * a[2] + sH * ml * Q * a[7] - H * N * a[4]);
		s[7][EQ] = d[7] - h * (E1 * a[3] - sH * mh * Q * a[6] - H * M * a[5]);
		s[8][EQ] = d[8] - h * (E2 * a[4] - sH * mh * P * a[9] - H * N * a[2]);
		s[9][EQ] = d[9] - h * (E3 * a[5] + sH * ml * P * a[8] - H * M * a[3]);

		s[10][EQ] = d[10];

		// Fill derivatives
		// Copy constant elements from pre-calculated D matrix first
		for (i = 1; i <= NE; ++i)
			for (j = 1; j <= NE; ++j)
				if (i != j)
					s[i][j] = s[i][NE+j] = h2 * D[i][j];

		// Fill solution-dependent elements 
		s[1][2] = s[1][NE+2] = -h * a[2];
		s[1][3] = s[1][NE+3] = -h * a[3];
		s[1][4] = s[1][NE+4] = -h * a[4];
		s[1][5] = s[1][NE+5] = -h * a[5];

		s[6][2] = s[6][NE+2] = -h2 * E0;
		s[7][3] = s[7][NE+3] = -h2 * E1;
		s[8][4] = s[8][NE+4] = -h2 * E2;
		s[9][5] = s[9][NE+5] = -h2 * E3;

		s[6][10] = s[6][NE+10] = h2 * a[2];
		s[7][10] = s[7][NE+10] = h2 * a[3];
		s[8][10] = s[8][NE+10] = h2 * a[4];
		s[9][10] = s[9][NE+10] = h2 * a[5];
	}
}

void SQWHSolver::init_derivative_matrix()
{
	clear_matrix(D, 1, NE, 1, NE);

	D[2][6] = -mh;
	D[3][7] = -ml;
	D[4][8] = -ml;
	D[5][9] = -mh;

	D[2][3] =  sH * mh * Q;
	D[3][2] = -sH * ml * Q;
	D[4][5] = -sH * ml * P;
	D[5][4] =  sH * mh * P;
	D[6][7] =  sH * ml * Q;
	D[7][6] = -sH * mh * Q;
	D[8][9] = -sH * mh * P;
	D[9][8] =  sH * ml * P;
	
	D[6][4] = -H * N;
	D[7][5] = -H * M;
	D[8][2] = -H * N;
	D[9][3] = -H * M;	
}

void SQWHSolver::eq_cb(int k, int* idx, float **s, float **y, void* ctx)
{
	static_cast<SQWHSolver*>(ctx)->eq(k, idx, s, y);
}

void SQWHSolver::solve_once(unsigned spin, float precision)
{
	solvde( eq_cb, ITMAX, precision, SLOWC, scalv, indexv, NE, NB, p.M, y, c, s, NULL, this );
	check( count_zeros(spin, precision) == p.subband, SQWH_QNCH );
}

void SQWHSolver::solve_zero_field()
{
	init_level(0);

	unsigned spin, i;
	for ( spin = 0; spin < 4; ++spin ) {
		init_guess(spin);
		for ( i = 1 ; i <= p.psteps ; i++ ) {
			pscale = float( i ) / p.psteps;
			solve_once( spin, CONV0 );
		}
		solve_once(spin, p.prec);
		copy_matrix(y, z[spin], 1, NE, 1, p.M, 1, 1);
	}
}

void SQWHSolver::solve_level(unsigned l)
{
	for ( unsigned spin = 0; spin < 4; ++spin ) {
		copy_matrix(z[spin], y, 1, NE, 1, p.M, 1, 1);
		init_level(l + spin);
		for ( unsigned b = 0; b <= p.Bsteps; ++b ) {
			set_magnetic_field(b * p.Bstep);
			solve_once(spin, p.prec);
			e[l][b][spin] = y[NE][1];
		}
	}
}

void SQWHSolver::save_wavefunction(float **f, const std::string& filename) const
{
	std::ofstream out( filename.c_str(), std::ios::trunc );
	check3( out, SQWH_OUT, filename );
	out << "\"E" << p.subband << '=' << f[NE][1] * p.Eo;
	for ( unsigned k = 1 ; k <= p.M ; k++ ) {
		out << std::endl
			<< ( k - 1.0 ) / ( p.M - 1.0 )  << ' '
			<< pot[k] * p.Eo << ' ' 
			<< f[2][k] << ' '
			<< f[3][k] << ' '
			<< f[4][k] << ' '
			<< f[5][k];
	}
	out.close();
}

void SQWHSolver::save_level(float **en, const std::string& filename) const
{
	std::ofstream out( filename.c_str(), std::ios::trunc );
	check3( out, SQWH_OUT, filename );
	for ( unsigned b = 0; b <= p.Bsteps; ++b ) {
		out << b * p.Bstep * p.Bo;
		for ( unsigned spin = 0; spin < 4; ++spin )
			out << ' ' << en[b][spin] * p.Eo;
		out << std::endl;
	}
}

