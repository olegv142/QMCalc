#include "sqwe.h"
#include "nrutil.h"
#include "nrtools.h"
#include "express.h"
#include "misc.h"
#include <fstream>  

#define _USE_MATH_DEFINES
#include <math.h>

#define ENE0 5
#define ENE1 4
#define ENB0 3
#define ENB1 2
#define NE ENE0
#define NC 3

#define PARSTART 25
#define TILTSTEP 3
#define PARSTEP .1f
#define ITMAX 1000
#define CONV0 .1f
#define SLOWC CONV0

DefineErrorMess( SQWE_OUT, "Can't write output file" );
DefineErrorMess( SQWE_QNCH, "Wavefunction quench" );
DefineErrorMess( SQWE_CONV, "Failed to converge" );

void GetSQWEParams(struct SQWEParams& p, TExpression const& e)
{
	p.M        = (unsigned)e.Get("grd");
	p.subbands = (unsigned)e.Get("subbands");
	p.psteps   = (unsigned)e.Get("psteps");
	p.W        = (float   )e.Get("W");
	p.me       = (float   )e.Get("me");
	p.me0      = (float   )e.Get("me0");
	p.me1      = (float   )e.Get("me1");
	p.ef       = (float   )e.Get("ef");
	p.e1       = (float   )e.Get("e1");
	p.cex      = (float   )e.Get("cex");
	p.eb0      = (float   )e.Get("eb0");
	p.eb1      = (float   )e.Get("eb1");
	p.Eo       = (float   )e.Get("Eo");
	p.prec     = (float   )e.Get("prec");
}

SQWESolver::SQWESolver(struct SQWEParams const& params)
	: p(params)
{
	assert( p.subbands > 0 && p.M > 1 && p.psteps > 1 );
	assert( p.me > 0 && p.me0 > 0 && p.me1 > 0 );
	assert( p.eb0 >= 0 && p.eb1 >= 0 );

	s = matrix( 1, NE, 1, 2 * NE + 1 );
	c = f3tensor( 1, NE, 1, NC, 1, p.M + 1 );
	yg = matrix( 1, ENE0, 1, p.M );
	if ( p.subbands > 1 )
		ys = f3tensor( 1, p.subbands - 1, 1, ENE1, 1, p.M );
	else 
		ys = NULL;

	pot = vector( 1, p.M );
	scalv0 = vector( 1, NE );
	scalvs = vector( 1, NE );
	indexv = ivector( 1, NE );

	unsigned i;
	for( i = 1; i <= NE ; i++ ) {
		scalv0[i] = scalvs[i] = 1;
		indexv[i] = i;
	}
	scalv0[ENE0] = scalvs[ENE1] = 1 / p.me;
	float B2 = ( 1 - p.W ) / 2;
	for( i = 1 ; i <= p.M ;  i++ ) {
		float z = ( i - 1.f ) / ( p.M - 1.f );
		if( z < B2 ) pot[i] = p.eb0;
		else if( z > 1 - B2 ) pot[i] = p.eb1;
		else pot[i] = 0;
	}
	pscale = 0;
}

SQWESolver::~SQWESolver()
{
	free_matrix( s, 1, NE, 1, 2 * NE + 1 );
	free_f3tensor( c, 1, NE, 1, NC, 1, p.M + 1 );
	free_matrix( yg, 1, ENE0, 1, p.M );
	if ( ys )
		free_f3tensor( ys, 1, p.subbands - 1, 1, ENE1, 1, p.M );
	free_vector( pot, 1, p.M );
	free_vector( scalv0, 1, NE );
	free_vector( scalvs, 1, NE );
	free_ivector( indexv, 1, NE );
}

// Build initial guess as sinusoidal function of z
void SQWESolver::eguess()
{
	assert( yg );
	assert( p.subbands <= 1 || ys );
	float const sq2 = sqrt( 2.f );
	float const pi  = (float)M_PI;
	unsigned k, i;
	for( k = 1 ; k <= p.M ; k++ ) {  
		float x = ( k - 1 ) / ( p.M - 1.0f );
		float px = pi * x;
		yg[2][k] = x - sin( 2 * px ) / ( 2 * pi );
		yg[1][k] = 0;
		yg[3][k] = sq2 * sin( px );
		yg[4][k] = pi * sq2 * cos( px ) / ( 2 * p.me );
		yg[5][k] = pi * pi / ( 2 * p.me );
		for( i = 1 ; i < p.subbands ; i++ ) {
			unsigned j = i + 1;
			float jpi = j * pi;
			float jpx = jpi * x;
			ys[i][1][k] = x - sin( 2 * jpx ) / ( 2 * jpi );
			ys[i][2][k] = sq2 * sin( jpx );
			ys[i][3][k] = jpi * sq2 * cos( jpx ) / ( 2 * p.me );
			ys[i][4][k] = jpi * jpi / ( 2 * p.me );
		}
	}
}

// Transform initial guess to rough approximation
void SQWESolver::eintro()
{
	assert( yg );
	assert( p.subbands <= 1 || ys );
	float eef = p.ef;
	float ee1 = p.e1;
	p.ef = p.e1 = 0;
	unsigned i;
	for( i = 1 ; i <= p.psteps ; i++ ) {
		pscale = float( i ) / p.psteps;
		esolve( CONV0 );
	}
	while( p.e1 < ee1 ) {  // apply electric field in steps
		esolve( CONV0 );
		p.e1 += TILTSTEP / p.me;
	}
	p.e1 = ee1;
	while( p.ef < eef ) {  // increase electron density in steps
		esolve( CONV0 );
		p.ef += PARSTART / p.me + PARSTEP * p.ef;
	}
	p.ef = eef;
}

void SQWESolver::setEnergy(unsigned subband, float val)
{
	assert( 0 <= subband && subband < p.subbands );
	float *e = ( subband > 0 ? ys[subband][4] : yg[5] );
	for( unsigned i = 1 ; i <= p.M ; i++ ) e[i] = val;
}

float SQWESolver::getEnergy(unsigned subband) const
{
	assert( 0 <= subband && subband < p.subbands );
	return subband == 0 ? yg[5][1] : ys[subband][4][1];
}

const float* SQWESolver::WaveFunction( unsigned subband ) const
{
	assert( 0 <= subband && subband < p.subbands );
	return subband == 0 ? &yg[3][1] : &ys[subband][2][1];
}

void SQWESolver::esolve(float precision)
{
	assert( yg );
	assert( p.subbands <= 1 || ys );
	if (solvde( eq0cb, ITMAX, precision, SLOWC, scalv0, indexv, ENE0, ENB0, p.M, yg, c, s, NULL, this ) > ITMAX)
		Signal( SQWE_CONV );
	checkZeroes( 0, precision );
	for( unsigned i = 1 ; i < p.subbands ; i++ ) {
		if (solvde( eq1cb, ITMAX, precision, SLOWC, scalvs, indexv, ENE1, ENB1, p.M, ys[i], c, s, NULL, this ) > ITMAX)
			Signal( SQWE_CONV );
		checkZeroes( i, precision );
	}
}

void SQWESolver::checkZeroes( unsigned subband, float precision ) const
{
	const float *f = WaveFunction( subband );
	unsigned cnt = count_zeros(f, 0, p.M - 1, precision);
	check( cnt == subband, SQWE_QNCH );
}

void SQWESolver::eq0cb(int k, int* idx, float **s, float **y, void* ctx)
{
	static_cast<SQWESolver*>(ctx)->eq0(k, idx, s, y);
}

void SQWESolver::eq1cb(int k, int* idx, float **s, float **y, void* ctx)
{
	static_cast<SQWESolver*>(ctx)->eq1(k, idx, s, y);
}

void SQWESolver::eq0(int k, int* idx, float **s, float **y) const
{
	float h = 1 / (p.M - 1.f);
	float h2 = h / 2;
	clear_matrix( s, 1, ENE0, 1, 2 * ENE0 + 1 );

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
		float a[ENE0+1];
		float d[ENE0+1];

		for( int i = 1 ; i <= ENE0 ; i++ ) {
			a[i] = ( y[i][k] + y[i][k-1] ) / 2;
			d[i] = y[i][k] - y[i][k-1];
			s[i][i] = -1;
			s[i][i+ENE0] = 1;
		}

		float a32 = a[3] * a[3];
		float pex = p.cex * pow( a32, 1.0f / 3 );
		float epot = pscale * ( pot[k] + pot[k-1] ) / 2;

		s[1][2] = s[1][7]  = p.ef * h2;
		s[2][3] = s[2][8]  = -h * a[3];
		s[3][4] = s[3][9]  = -h * p.me;
		s[4][1] = s[4][6]  = -h2 * a[3];
		s[4][3] = s[4][8]  = h2 * ( a[5] - a[1] - epot  + ( 5.0f / 3 ) * pex );
		s[4][5] = s[4][10] = h2 * a[3];

		s[1][11] = d[1] + h * ( p.ef * ( a[2] - 1 ) - p.e1 );
		s[2][11] = d[2] - h * a32;
		s[3][11] = d[3] - 2 * h * a[4] * p.me;
		s[4][11] = d[4] + h * ( a[5] - a[1] - epot + pex ) * a[3];
		s[5][11] = d[5];
	}
}

void SQWESolver::eq1(int k, int* idx, float **s, float **y) const
{
	float h = 1 / (p.M - 1.f);
	float h2 = h / 2;
	clear_matrix( s, 1, ENE1, 1, 2 * ENE1 + 1 );

	if( k == 1 ) {
		s[3][5] = s[4][6] = 1;
		s[3][9] = y[1][1];
		s[4][9] = y[2][1];
	} else if( k > (int)p.M ) {
		s[1][5] = s[2][6] = 1;
		s[1][9] = y[1][p.M] - 1;
		s[2][9] = y[2][p.M];
	} else {
		float a[ENE1+1];
		float d[ENE1+1];

		for( int i = 1 ; i <= ENE1 ; i++ ) {
			a[i] = ( y[i][k] + y[i][k-1] ) / 2;
			d[i] = y[i][k] - y[i][k-1];
			s[i][i] = -1;
			s[i][i+ENE1] = 1;
		}

		float epot = pscale * ( pot[k] + pot[k-1] ) / 2;
		float const * const elpot = yg[1];
		a[0] = ( elpot[k] + elpot[k-1] ) / 2;

		s[1][2] = s[1][6] = -h * a[2];
		s[2][3] = s[2][7] = -h * p.me;
		s[3][2] = s[3][6] = h2 * ( a[4] - a[0] - epot );
		s[3][4] = s[3][8] = h2 * a[2];

		s[1][9] = d[1] - h * a[2] * a[2];
		s[2][9] = d[2] - 2 * h * a[3] * p.me;
		s[3][9] = d[3] + h * ( a[4] - a[0] - epot ) * a[2];
		s[4][9] = d[4];
	}
}

void SQWESolver::Solve()
{
	eguess(); // Make initial guess
	eintro(); // Tune it to account electric field
	esolve( p.prec ); // Solve
}

void SQWESolver::SaveResults(const char* filename) const
{
	std::ofstream out( filename, std::ios::trunc );
	check3( out, SQWE_OUT, filename );
	out << "\"E0=" << getEnergy( 0 ) * p.Eo << "; ";
	for( unsigned i = 1 ; i < p.subbands ; i++ ) 
		out << "E" << i << "0=" << ( getEnergy( i ) - getEnergy( 0 ) ) * p.Eo << "; ";
	for( unsigned k = 1 ; k <= p.M ; k++ ) {
		out << std::endl;
		out << ( k - 1.0 ) / ( p.M - 1.0 )  << " " 
			<<  ( yg[1][k] + pot[k] ) * p.Eo << " " 
			<< yg[3][k];
		for( unsigned i = 1 ; i < p.subbands ; i++ ) out << " " <<  ys[i][2][k];
	}
	out.close();
}
