#include "sqwh.h"
#include "nrutil.h"
#include "nrtools.h"
#include "express.h"
#include "misc.h"
#include <fstream>  

#define _USE_MATH_DEFINES
#include <math.h>

#define ITMAX 512
#define SLOWC 10.f

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
DefineErrorMess( SQWH_CONV, "Failed to converge" );

void GetSQWHParams(struct SQWHParams& p, TExpression const& e)
{
	p.M        = (unsigned)e.Get("grid");
	p.NL       = (unsigned)e.Get("NL");
	p.subband  = (unsigned)e.Get("subband");
	p.honly    = (0 != e.Get("honly"));
	p.hb       = (float   )e.Get("hb");
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

	B0 = matrix( 1, 8, 1, 8 );
	B1 = matrix( 1, 8, 1, 8 );

	for (unsigned i = 0; i < 4; ++i) {
		sol_z[i] = matrix( 1, NE, 1, p.M );
		sol_h[i] = f3tensor( 0, p.NL - 1, 1, NE, 1, p.M );
	}

	if (p.NL > 0) {
		sol_e = f3tensor( 0, p.NL - 1, 0, p.Bsteps, 0, 3 );
		clear_f3tensor( sol_e, 0, p.NL - 1, 0, p.Bsteps, 0, 3 );
	} else
		sol_e = NULL;
}

SQWHSolver::~SQWHSolver()
{
	free_matrix( y, 1, NE, 1, p.M );
	free_matrix( s, 1, NE, 1, 2 * NE + 1 );
	free_f3tensor( c, 1, NE, 1, NE - NB + 1, 1, p.M + 1 );
	free_vector( scalv, 1, NE );
	free_ivector( indexv, 1, NE );
	free_matrix( B0, 1, 8, 1, 8 );
	free_matrix( B1, 1, 8, 1, 8 );
	for (unsigned i = 0; i < 4; ++i) {
		free_matrix( sol_z[i], 1, NE, 1, p.M );
		free_f3tensor( sol_h[i], 0, p.NL - 1, 1, NE, 1, p.M );
	}
	if (p.NL > 0)
		free_f3tensor( sol_e, 0, p.NL - 1, 0, p.Bsteps, 0, 3 );
}

// Create initial guess for the particular spin s = 0..3
void SQWHSolver::init_guess(unsigned spin)
{
	float m = (spin == 0 || spin == 3) ? mh : ml;
	float pl = (float)M_PI * (p.subband + 1.f);
	float e = pl * pl / m;
	clear_matrix( y, 1, NE, 1, p.M );
	for( unsigned i = 1 ; i <= p.M ;  i++ ) {
		float z = ( i - 1.f ) / ( p.M - 1.f );
		float f = pl * z;
		y[1][i] = z - sin(2 * f) / (2 * pl);
		y[2 + spin][i] = sqrt(2.f) * sin(f);
		y[6 + spin][i] = sqrt(2.f) * pl * cos(f) / m;
		y[10][i] = e;
	}
}

void SQWHSolver::set_params(unsigned n, float field)
{
	float n0 = n + .5f;
	float n1 = n > 0 ? n - 1 + .5f : 0;
	float n2 = n > 1 ? n - 2 + .5f : 0;
	float n3 = n > 2 ? n - 3 + .5f : 0;
	float R = sqrt(3.f) * (p.g2 + p.g3) / 2;
	float S = sqrt(3/2.f) * p.g3;

	M = n > 2 ? R * sqrt((n - 1.f)*(n - 2.f)) : 0;
	N = n > 1 ? R * sqrt(n*(n - 1.f)) : 0;
	P = n > 2 ? S * sqrt(n - 2.f) : 0;
	Q = n > 0 ? S * sqrt((float)n) : 0;

	H = field;
	sH = sqrt(field);

	E0 = (A*n0 - ml*Q*Q - 3*p.K/2)*H;
	E1 = (B*n1 - mh*Q*Q -   p.K/2)*H;
	E2 = (B*n2 - mh*P*P +   p.K/2)*H;
	E3 = (A*n3 - ml*P*P + 3*p.K/2)*H;
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

	save_wavefunction(sol_z[0], bname + "Z_hm.dat");
	if (!p.honly) {
		save_wavefunction(sol_z[1], bname + "Z_lm.dat");
		save_wavefunction(sol_z[2], bname + "Z_lp.dat");
	}
	save_wavefunction(sol_z[3], bname + "Z_hp.dat");

	bname += 'L';
	for (unsigned l = 0; l < p.NL; ++l) {
		save_wavefunction(sol_h[0][l], bname + (char)('0' + l) + "_hm.dat");
		if (!p.honly) {
			save_wavefunction(sol_h[1][l], bname + (char)('0' + l) + "_lm.dat");
			save_wavefunction(sol_h[2][l], bname + (char)('0' + l) + "_lp.dat");
		}
		save_wavefunction(sol_h[3][l], bname + (char)('0' + l) + "_hp.dat");
	}

	save_levels(bname + ".dat");
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
	int i;
	clear_matrix( s, 1, NE, 1, 2 * NE + 1 );
#define EQ (2*NE+1)
	if( k == 1 ) {
		// The y[1] = 0 at left boundary
		s[6][EQ] = y[1][1];
		s[6][NE+1] = 1;
		for (i = 1; i <= 8; ++i) {
			float f = y[1+i][1];
			s[7][EQ]  += f * B0[1][i];
			s[8][EQ]  += f * B0[2][i];
			s[9][EQ]  += f * B0[3][i];
			s[10][EQ] += f * B0[4][i];
			s[7][NE+1+i]  = B0[1][i];
			s[8][NE+1+i]  = B0[2][i];
			s[9][NE+1+i]  = B0[3][i];
			s[10][NE+1+i] = B0[4][i];
		}
	} else if( k > (int)p.M ) {
		// The y[1] = 1 at the right boundary
		s[1][EQ]  = y[1][p.M] - 1;
		s[1][NE+1] = 1;
		for (i = 1; i <= 8; ++i) {
			float f = y[1+i][p.M];
			s[2][EQ] += f * B1[1][i];
			s[3][EQ] += f * B1[2][i];
			s[4][EQ] += f * B1[3][i];
			s[5][EQ] += f * B1[4][i];
			s[2][NE+1+i] = B1[1][i];
			s[3][NE+1+i] = B1[2][i];
			s[4][NE+1+i] = B1[3][i];
			s[5][NE+1+i] = B1[4][i];
		}
	} else {
		float h = 1 / (p.M - 1.f), h2 = h / 2;
		// Internal pair of points
		float a[NE+1]; // 2 point average
		float d[NE+1]; // 2 point deltas

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
		s[6][EQ] = d[6] + h * ((a[10]-E0) * a[2] + sH * ml * Q * a[7] - H * N * a[4]);
		s[7][EQ] = d[7] + h * ((a[10]-E1) * a[3] - sH * mh * Q * a[6] - H * M * a[5]);
		s[8][EQ] = d[8] + h * ((a[10]-E2) * a[4] - sH * mh * P * a[9] - H * N * a[2]);
		s[9][EQ] = d[9] + h * ((a[10]-E3) * a[5] + sH * ml * P * a[8] - H * M * a[3]);
		s[10][EQ] = d[10];

		// Fill derivatives
		get_derivatives(a, s, 0,  h2);
		get_derivatives(a, s, NE, h2);
	}
}

void SQWHSolver::get_derivatives(float *a, float **D, int shift, float mult) const
{
	D[1][shift+2] = -2 * a[2] * mult;
	D[1][shift+3] = -2 * a[3] * mult;
	D[1][shift+4] = -2 * a[4] * mult;
	D[1][shift+5] = -2 * a[5] * mult;

	D[2][shift+6] = -mh * mult;
	D[3][shift+7] = -ml * mult;
	D[4][shift+8] = -ml * mult;
	D[5][shift+9] = -mh * mult;

	D[2][shift+3] =  sH * mh * Q * mult;
	D[3][shift+2] = -sH * ml * Q * mult;
	D[4][shift+5] = -sH * ml * P * mult;
	D[5][shift+4] =  sH * mh * P * mult;
	D[6][shift+7] =  sH * ml * Q * mult;
	D[7][shift+6] = -sH * mh * Q * mult;
	D[8][shift+9] = -sH * mh * P * mult;
	D[9][shift+8] =  sH * ml * P * mult;
	
	D[6][shift+4] = -H * N * mult;
	D[7][shift+5] = -H * M * mult;
	D[8][shift+2] = -H * N * mult;
	D[9][shift+3] = -H * M * mult;

	D[6][shift+2] = (a[10]-E0) * mult;
	D[7][shift+3] = (a[10]-E1) * mult;
	D[8][shift+4] = (a[10]-E2) * mult;
	D[9][shift+5] = (a[10]-E3) * mult;

	D[6][shift+10] = a[2] * mult;
	D[7][shift+10] = a[3] * mult;
	D[8][shift+10] = a[4] * mult;
	D[9][shift+10] = a[5] * mult;
}

// Get equation matrix for the particular potential value
void SQWHSolver::get_equation(float **m, float pot)
{
	clear_matrix(m, 1, 8, 1, 8);

	// Build matrix for equation y' + m*y = 0
	m[1][5] = m[4][8] = -mh;
	m[2][6] = m[3][7] = -ml;
	m[1][2] =  sH*mh*Q;
	m[2][1] = -sH*ml*Q;
	m[3][4] = -sH*ml*P;
	m[4][3] =  sH*mh*P;
	m[5][6] =  sH*ml*Q;
	m[6][5] = -sH*mh*Q;
	m[7][8] = -sH*mh*P;
	m[8][7] =  sH*ml*P;
	m[5][3] = -H*N;
	m[6][4] = -H*M;
	m[7][1] = -H*N;
	m[8][2] = -H*M;

	float E = y[NE][1] - pot;
	m[5][1] = E - E0;
	m[6][2] = E - E1;
	m[7][3] = E - E2;
	m[8][4] = E - E3;
}

void SQWHSolver::init_boundary_condition()
{
	int i, j;
	float **m = matrix(1, 8, 1, 8);
	get_equation(m, p.hb);
	emod(m, 8, p.prec, ITMAX, B0);
	copy_matrix(B0, B1, 1, 8, 1, 8, 1, 1);
	for (i = 1; i <= 8; ++i)
		for (j = 1; j <= 8; ++j) {
			B0[i][j] += m[i][j];
			B1[i][j] -= m[i][j];
		}

	free_matrix(m, 1, 8, 1, 8);
}

void SQWHSolver::eq_cb(int k, int* idx, float **s, float **y, void* ctx)
{
	static_cast<SQWHSolver*>(ctx)->eq(k, idx, s, y);
}

void SQWHSolver::solve_once(unsigned spin, float precision)
{
	for (int it_max = 1;; it_max *= 2) {
		init_boundary_condition();
		if (solvde( eq_cb, it_max, precision, SLOWC * p.M, scalv, indexv, NE, NB, p.M, y, c, s, NULL, this ) <= it_max)
			break;
		if (it_max >= ITMAX)
			Signal(SQWH_CONV);
	}
}

void SQWHSolver::solve_zero_field()
{
	set_params(0, 0);
	for ( unsigned spin = 0; spin < 4; ++spin ) {
		if (skip_light_hole(spin))
			continue;
		init_guess(spin);
		solve_once(spin, p.prec);
		copy_matrix(y, sol_z[spin], 1, NE, 1, p.M, 1, 1);
	}
}

void SQWHSolver::solve_level(unsigned l)
{
	for ( unsigned spin = 0; spin < 4; ++spin ) {
		if (skip_light_hole(spin))
			continue;
		copy_matrix(sol_z[spin], y, 1, NE, 1, p.M, 1, 1);
		for ( unsigned b = 0; b <= p.Bsteps; ++b ) {
			set_params(l + spin, b * p.Bstep);
			if (b > 1) {
				float E = 2*sol_e[l][b-1][spin] - sol_e[l][b-2][spin];
				for (unsigned i = 1; i <= p.M; ++i)
					y[NE][i] = E;
			}
			solve_once(spin, p.prec);
			sol_e[l][b][spin] = y[NE][1];
		}
		copy_matrix(y, sol_h[spin][l], 1, NE, 1, p.M, 1, 1);
	}
}

void SQWHSolver::save_wavefunction(float **f, const std::string& filename) const
{
	std::ofstream out( filename.c_str(), std::ios::trunc );
	check3( out, SQWH_OUT, filename );
	out << "\"E" << f[NE][1] * p.Eo;
	for ( unsigned k = 1 ; k <= p.M ; k++ ) {
		out << std::endl
			<< ( k - 1.0 ) / ( p.M - 1.0 )  << ' '
			<< f[2][k] << ' '
			<< f[3][k] << ' '
			<< f[4][k] << ' '
			<< f[5][k];
	}
	out.close();
}

void SQWHSolver::save_levels(const std::string& filename) const
{
	std::ofstream out( filename.c_str(), std::ios::trunc );
	check3( out, SQWH_OUT, filename );
	for ( unsigned b = 0; b <= p.Bsteps; ++b ) {
		out << b * p.Bstep * p.Bo;
		for (unsigned l = 0; l < p.NL; ++l)
			for ( unsigned spin = 0; spin < 4; ++spin ) {
				if (skip_light_hole(spin))
					continue;
				out << ' ' << sol_e[l][b][spin] * p.Eo;
			}
		out << std::endl;
	}
}

bool SQWHSolver::skip_light_hole(unsigned spin) const
{
	return p.honly && (spin == 1 || spin == 2);
}
