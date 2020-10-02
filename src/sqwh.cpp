#include "sqwh.h"
#include "nrutil.h"
#include "nrtools.h"
#include "express.h"
#include "misc.h"
#include <fstream>  

#define _USE_MATH_DEFINES
#include <math.h>

#define ITMAX 512
#define CONV0 .1f
#define SLOWC CONV0

// The number of wavefunction components
#define NF 8
// The number of equations
#define NE (NF+2)
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
	p.psteps   = (unsigned)e.Get("psteps");
	p.hb       = (float   )e.Get("hb");
	p.Bstep    = (float   )e.Get("Bstep");
	p.Bsteps   = (unsigned)e.Get("Bsteps");
	p.Bskip    = (unsigned)e.Get("Bskip");
	p.ef       = (float   )e.Get("ef");
	p.me       = (float   )e.Get("me");
	p.Eo       = (float   )e.Get("Eo");
	p.Bo       = (float   )e.Get("Bo");
	p.g1       = (float   )e.Get("g1");
	p.g2       = (float   )e.Get("g2");
	p.g3       = (float   )e.Get("g3");
	p.K        = (float   )e.Get("K");
	p.e_side   = (float   )e.Get("e_side");
	p.e_cyc    = (float   )e.Get("e_cyc");
	p.e_spin   = (float   )e.Get("e_spin");
	p.prec     = (float   )e.Get("prec");
}

SQWHSolver::SQWHSolver(struct SQWHParams const& params)
	: p(params)
{
	assert(p.NL > 0);
	assert(p.Bskip > 0);

	y = matrix( 1, NE, 1, p.M );
	s = matrix( 1, NE, 1, 2 * NE + 1 );
	c = f3tensor( 1, NE, 1, NE - NB + 1, 1, p.M + 1 );
	pot = vector( 1, p.M );
	scalv  = vector( 1, NE );
	indexv = ivector( 1, NE );

	mh = 1 / (p.g1 / 2 - p.g2);
	ml = 1 / (p.g1 / 2 + p.g2);
	A = p.g1 + p.g2;
	B = p.g1 - p.g2;

	B0 = matrix( 1, NF, 1, NF );
	B1 = matrix( 1, NF, 1, NF );

	sol_el = vector( 1, p.M );
	for (unsigned i = 0; i < 4; ++i) {
		sol_z[i] = matrix( 1, NE, 1, p.M );
		sol_f[i] = f4tensor( 0, p.NL - 1, 0, p.Bsteps, 1, NF, 1, p.M );
		sol_e[i] = matrix( 0, p.NL - 1, 0, p.Bsteps);
		clear_matrix( sol_e[i], 0, p.NL - 1, 0, p.Bsteps);
	}
}

SQWHSolver::~SQWHSolver()
{
	free_matrix( y, 1, NE, 1, p.M );
	free_matrix( s, 1, NE, 1, 2 * NE + 1 );
	free_f3tensor( c, 1, NE, 1, NE - NB + 1, 1, p.M + 1 );
	free_vector( pot, 1, p.M );
	free_vector( scalv, 1, NE );
	free_ivector( indexv, 1, NE );
	free_matrix( B0, 1, NF, 1, NF );
	free_matrix( B1, 1, NF, 1, NF );
	free_vector( sol_el, 1, p.M );
	for (unsigned i = 0; i < 4; ++i) {
		free_matrix( sol_z[i], 1, NE, 1, p.M );
		free_f4tensor( sol_f[i], 0, p.NL - 1, 0, p.Bsteps, 1, NF, 1, p.M );
		free_matrix( sol_e[i], 0, p.NL - 1, 0, p.Bsteps);
	}
}

#define ENE 5
#define ENB 3
// Electron equations:
// y[1][1..M] - electrostatic potential
// y[2][1..M] - wavefunction density (psy^2) integral
// y[3][1..M] - wavefunction
// y[4][1..M] - wavefunction derivative
// y[5][1..M] - energy eigenvalue

void SQWHSolver::init_eguess()
{
	float const sq2 = sqrt( 2.f );
	float const pi  = (float)M_PI;
	unsigned i;
	for( i = 1; i <= NE ; i++ ) {
		scalv[i] = 1;
		indexv[i] = i;
	}
	scalv[ENE] = pi * pi / ( 2 * p.me );
	for( i = 1 ; i <= p.M ; i++ ) {  
		float x = ( i - 1 ) / ( p.M - 1.0f );
		float px = pi * x;
		y[1][i] = 0;
		y[2][i] = x - sin( 2 * px ) / ( 2 * pi );
		y[3][i] = sq2 * sin( px );
		y[4][i] = pi * sq2 * cos( px ) / ( 2 * p.me );
		y[5][i] = pi * pi / ( 2 * p.me );
	}
}

void SQWHSolver::eqe(int k, int* idx, float **s, float **y) const
{
	float h = 1 / (p.M - 1.f);
	float h2 = h / 2;
	float ef = p.ef * pscale;
	clear_matrix( s, 1, ENE, 1, 2 * ENE + 1 );

	if( k == 1 ) {
		s[3][6] = s[4][7] = s[5][8] = 1;
		s[3][11] = y[1][1];
		s[4][11] = y[2][1];
		s[5][11] = y[3][1];
	} else if( k > (int)p.M ) {
		s[1][7] = s[2][8] = 1;
		s[1][11] = y[2][p.M] - 1;
		s[2][11] = y[3][p.M];
	} else {
		float a[ENE+1];
		float d[ENE+1];

		for( int i = 1 ; i <= ENE ; i++ ) {
			a[i] = ( y[i][k] + y[i][k-1] ) / 2;
			d[i] = y[i][k] - y[i][k-1];
			s[i][i] = -1;
			s[i][i+ENE] = 1;
		}

		s[1][2] = s[1][7]  = ef * h2;
		s[2][3] = s[2][8]  = -h * a[3];
		s[3][4] = s[3][9]  = -h * p.me;
		s[4][1] = s[4][6]  = -h2 * a[3];
		s[4][3] = s[4][8]  = h2 * ( a[5] - a[1] );
		s[4][5] = s[4][10] = h2 * a[3];

		s[1][11] = d[1] + h * ef * ( a[2] - 1 + p.e_side );
		s[2][11] = d[2] - h * a[3] * a[3];
		s[3][11] = d[3] - 2 * h * a[4] * p.me;
		s[4][11] = d[4] + h * ( a[5] - a[1] ) * a[3];
		s[5][11] = d[5];
	}
}

void SQWHSolver::solve_electrons()
{
	init_eguess();
	for (unsigned step = 0; step <= p.psteps; ++step) {
		pscale = (float)step / p.psteps;
		solve_e(CONV0);
	}
	solve_e(p.prec);
	copy_vector(y[1], pot, 1, p.M, 1);
	copy_vector(y[3], sol_el, 1, p.M, 1);
}

// Create initial guess for the particular spin s = 0..3
void SQWHSolver::init_guess(unsigned spin)
{
	float m = (spin == 0 || spin == 3) ? mh : ml;
	float pl = (float)M_PI * (p.subband + 1.f);
	float e = pl * pl / m;
	unsigned i;
	for( i = 1; i <= NE ; i++ ) {
		scalv[i] = 1;
		indexv[i] = i;
	}
	scalv[NE] = pl * pl / (ml + mh);
	clear_matrix( y, 1, NE, 1, p.M );
	for( i = 1 ; i <= p.M ;  i++ ) {
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
	solve_electrons();
	solve_zero_field();
	for (unsigned l = 0; l < p.NL; ++l)
		solve_level(l);
}

void SQWHSolver::SaveResults(const char* basename) const
{
	std::string bname(basename);
	bname += (char)('0' + p.subband);
	save_wavefunctions(0, 0, bname + 'Z');
	bname += 'L';
	for (unsigned l = 0; l < p.NL; ++l)
		save_wavefunctions(l, p.Bsteps, bname + (char)('0' + l));
	save_levels(bname + ".dat", false);
	save_levels(bname + "x.dat", true);
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
void SQWHSolver::eqh(int k, int* idx, float **s, float **y) const
{
	int i;
	clear_matrix( s, 1, NE, 1, 2 * NE + 1 );
#define EQ (2*NE+1)
	if( k == 1 ) {
		// The y[1] = 0 at left boundary
		s[6][EQ] = y[1][1];
		s[6][NE+1] = 1;
		for (i = 1; i <= NF; ++i) {
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
		for (i = 1; i <= NF; ++i) {
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
		// The electrostatic potential. It was positive for electrons.
		float hpot = -pscale * (pot[k] + pot[k-1]) / 2;
		// Internal pair of points
		float a[1+NE]; // 2 point average
		float d[1+NE]; // 2 point deltas

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
		s[6][EQ] = d[6] + h * ((a[10]-E0-hpot) * a[2] + sH * ml * Q * a[7] - H * N * a[4]);
		s[7][EQ] = d[7] + h * ((a[10]-E1-hpot) * a[3] - sH * mh * Q * a[6] - H * M * a[5]);
		s[8][EQ] = d[8] + h * ((a[10]-E2-hpot) * a[4] - sH * mh * P * a[9] - H * N * a[2]);
		s[9][EQ] = d[9] + h * ((a[10]-E3-hpot) * a[5] + sH * ml * P * a[8] - H * M * a[3]);
		s[10][EQ] = d[10];

		// Fill derivatives
		get_derivatives(a, hpot, s, 0,  h2);
		get_derivatives(a, hpot, s, NE, h2);
	}
}

void SQWHSolver::get_derivatives(float *a, float hpot, float **D, int shift, float mult) const
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

	D[6][shift+2] = (a[10]-E0-hpot) * mult;
	D[7][shift+3] = (a[10]-E1-hpot) * mult;
	D[8][shift+4] = (a[10]-E2-hpot) * mult;
	D[9][shift+5] = (a[10]-E3-hpot) * mult;

	D[6][shift+10] = a[2] * mult;
	D[7][shift+10] = a[3] * mult;
	D[8][shift+10] = a[4] * mult;
	D[9][shift+10] = a[5] * mult;
}

// Get equation matrix for the particular energy value
void SQWHSolver::get_equation(float **m, float E) const
{
	clear_matrix(m, 1, NF, 1, NF);

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
	m[5][1] = E - E0;
	m[6][2] = E - E1;
	m[7][3] = E - E2;
	m[8][4] = E - E3;
}

void SQWHSolver::init_boundary_condition()
{
	int i, j;
	float **m = matrix(1, NF, 1, NF);

	get_equation(m, y[NE][1] - p.hb + pscale * pot[1]);
	emod(m, NF, p.prec, ITMAX, B0);
	for (i = 1; i <= NF; ++i)
		for (j = 1; j <= NF; ++j)
			B0[i][j] += m[i][j];

	get_equation(m, y[NE][p.M] - p.hb + pscale * pot[p.M]);
	emod(m, NF, p.prec, ITMAX, B1);
	for (i = 1; i <= NF; ++i)
		for (j = 1; j <= NF; ++j)
			B1[i][j] -= m[i][j];

	free_matrix(m, 1, NF, 1, NF);
}

float SQWHSolver::barrier_integral(float const fb[NF], float E) const
{
	float I = 0;
	float **m = matrix(1, NF, 1, NF);
	float **mp = matrix(1, NF, 1, NF);
	float f[NF];
	float h, max = 0;
	int i, j;
	for (i = 0; i < NF; ++i)
		f[i] = fb[i];

	get_equation(m, E);
	// Get rid of the negative eigenvalues to ensure the solution decay
	emod(m, NF, p.prec, ITMAX, mp);
	// Choose step based on the maximum matrix element
	for (i = 1; i <= NF; ++i)
		for (j = 1; j <= NF; ++j)
			if (mp[i][j] > max)
				max = mp[i][j];

	h = sqrt(p.prec) / max;

	// Integrate wave function modulo squared while stepping into the barrier
	for (;;) {
		float f_[NF];
		for (i = 0; i < NF; ++i)
			f_[i] = f[i];
		for (i = 0; i < NF; ++i)
			for (j = 0; j < NF; ++j)
				f[i] -= mp[1+i][1+j] * f_[j] * h;
		float f2 = f[0] * f[0] + f[1] * f[1] + f[2] * f[2] + f[3] * f[3];
		I += f2 * h;
		if (f2 < p.prec)
			break;
	}

	free_matrix(m, 1, NF, 1, NF);
	free_matrix(mp, 1, NF, 1, NF);
	return I;
}

void SQWHSolver::solution_normalize(float **f, float e)
{
	float fb[NF], I = 0;
	unsigned i, j, side;
	for (side = 0; side < 2; ++side) {
		int k = side ? p.M : 1;
		float hpot = -pscale * pot[k];
		for (i = 0; i < NF; ++i)
			fb[i] = f[1+i][k];
		I += barrier_integral(fb, e - p.hb - hpot);
	}
	float mult = 1.f / (1.f + I);
	for (i = 1; i <= NF; ++i)
		for (j = 1; j <= p.M; ++j)
				f[i][j] *= mult;
}

void SQWHSolver::eqe_cb(int k, int* idx, float **s, float **y, void* ctx)
{
	static_cast<SQWHSolver*>(ctx)->eqe(k, idx, s, y);
}

void SQWHSolver::eqh_cb(int k, int* idx, float **s, float **y, void* ctx)
{
	static_cast<SQWHSolver*>(ctx)->eqh(k, idx, s, y);
}

void SQWHSolver::solve_e(float precision)
{
	if (solvde( eqe_cb, ITMAX, precision, SLOWC, scalv, indexv, ENE, ENB, p.M, y, c, s, NULL, this ) >= ITMAX)
		Signal(SQWH_CONV);
	unsigned cnt = count_zeros(y[3], 1, p.M, precision);
	check( cnt == 0, SQWH_QNCH );
}

void SQWHSolver::solve_h(unsigned spin, float precision)
{
	for (int it_max = 1;; it_max *= 2) {
		init_boundary_condition();
		if (solvde( eqh_cb, it_max, precision, SLOWC, scalv, indexv, NE, NB, p.M, y, c, s, NULL, this ) <= it_max)
			break;
		if (it_max >= ITMAX)
			Signal(SQWH_CONV);
	}
}

void SQWHSolver::solve_zero_field()
{
	set_params(0, 0);
	for ( unsigned spin = 0; spin < 4; ++spin ) {
		init_guess(spin);
		for (unsigned step = 0; step <= p.psteps; ++step) {
			pscale = (float)step / p.psteps;
			solve_h(spin, CONV0);
		}
		solve_h(spin, p.prec);
		copy_matrix(y, sol_z[spin], 1, NE, 1, p.M, 1, 1);
	}
}

void SQWHSolver::solve_level(unsigned l)
{
	for ( unsigned spin = 0; spin < 4; ++spin ) {
		copy_matrix(sol_z[spin], y, 1, NE, 1, p.M, 1, 1);
		for ( unsigned b = 0; b <= p.Bsteps; ++b ) {
			set_params(l + spin, b * p.Bstep);
			if (b > 1) {
				float E = 2*sol_e[spin][l][b-1] - sol_e[spin][l][b-2];
				for (unsigned i = 1; i <= p.M; ++i)
					y[NE][i] = E;
			}
			solve_h(spin, p.prec);
			sol_e[spin][l][b] = y[NE][1];
			copy_matrix(y + 1, sol_f[spin][l][b], 1, NF, 1, p.M, 1, 1);
			solution_normalize(sol_f[spin][l][b], y[NE][1]);
			printf("*");
		}
	}
}

void SQWHSolver::save_wavefunctions(unsigned n, unsigned b, const std::string& bname) const
{
	save_wavefunction(sol_f[0][n][b], bname + "_hm.dat");
	save_wavefunction(sol_f[1][n][b], bname + "_lm.dat");
	save_wavefunction(sol_f[2][n][b], bname + "_lp.dat");
	save_wavefunction(sol_f[3][n][b], bname + "_hp.dat");
}

void SQWHSolver::save_wavefunction(float **f, const std::string& filename) const
{
	std::ofstream out( filename.c_str(), std::ios::trunc );
	check3( out, SQWH_OUT, filename );
	for ( unsigned k = 1 ; k <= p.M ; k++ ) {
		out << ( k - 1.0 ) / ( p.M - 1.0 )  << ' ' << -pot[k] * p.Eo << ' '
			<< f[1][k] << ' '
			<< f[2][k] << ' '
			<< f[3][k] << ' '
			<< f[4][k] << std::endl;
	}
	out.close();
}

float SQWHSolver::get_intensity(float **f, unsigned n0) const
{
	if (n0 > 3)
		return 0;
	float I = 0;
	for( unsigned k = 1 ; k <= p.M ;  k++ )
		I += sol_el[k] * f[1+n0][k];
	if (n0 == 0 || n0 == 3)
		I *= 3;
	// Normalize so that the maximum intensity (heavy hole) will be 1;
	return fabs(I) / (3 * p.M);
}

void SQWHSolver::save_levels(const std::string& filename, bool transitions) const
{
	const char* spin_label[] = {"hm", "lm", "lp", "hp"};
	std::ofstream out( filename.c_str(), std::ios::trunc );
	check3( out, SQWH_OUT, filename );
	out << "\"B[T]";
	for (unsigned l = 0; l < p.NL; ++l)
		for ( unsigned spin = 0; spin < 4; ++spin ) {
			out << " L" << (char)('0' + l) << spin_label[spin];
			if (transitions) {
				out << " IL" << (char)('0' + l) << spin_label[spin]
					<< ((l + spin) % 4 < 2 ? "s-" : "s+") // polarization
					<< ((l + spin) % 2 < 1 ? "e-" : "e+");// electron spin
			}
		}
	out << std::endl;
	for ( unsigned b = 0; b <= p.Bsteps; b += p.Bskip ) {
		out << b * p.Bstep * p.Bo;
		for (unsigned l = 0; l < p.NL; ++l)
			for ( unsigned spin = 0; spin < 4; ++spin ) {
				float E = sol_e[spin][l][b];
				if (transitions)
					E += b * p.Bstep * (p.e_cyc + (((l + spin) % 2) - .5f) * p.e_spin);
				out << ' ' << E * p.Eo;
				if (transitions)
					out << ' ' << get_intensity(sol_f[spin][l][b], l + spin);
			}
		out << std::endl;
	}
}
