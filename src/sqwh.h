#pragma once

#include "errs.h"

// Parameters definitions for single quantum well holes model (sqwh.cpp)
// Units are selected to be dimensionless, see also comments in sqwh.par
struct SQWHParams {
	unsigned M;          // The number of grid points (in z direction)
	unsigned NL;         // Number of cyclotron levels computed
	unsigned subband;    // The quantized subband computed (zero based)
	unsigned psteps;     // Number of steps for potential turning on
	float    W;          // SQW width in the units of model region width
	float    hb0;        // Barriers height to the left  of the well (in dimensionless units)
	float    hb1;        // Barriers height to the right of the well (in dimensionless units)
	float    Bstep;      // Magnetic field step (in dimensionless units)
	unsigned Bsteps;     // Magnetic field steps (so B varies from 0 to B_step * B_steps)
	float    Eo;         // Energy scale in output files
	float    Bo;         // Magnetic field scale in output files
	float g1, g2, g3, K; // Luttinger parameters
};

class TExpression;

void GetSQWHParams(struct SQWHParams& p, TExpression const& e);

ErrorCode( SQWH_OUT );
ErrorCode( SQWH_QNCH );

class SQWHSolver {
public:
	SQWHSolver(struct SQWHParams const& params);
	~SQWHSolver();

	void  Solve(float precision);
	void  SaveResults(const char* basename) const;

protected:
	void init_guess(unsigned s);
	void init_level(unsigned lvl);
	void set_magnetic_field(float field);
	void init_derivative_matrix();
	unsigned count_zeros(unsigned s) const;
	void solve_once(unsigned spin, float precision);
	void solve_zero_field(float precision);
	void save_wavefunction(float **f, const char* filename) const;

	void eq(int k, int* idx, float **s, float **y) const;
	static void eq_cb(int k, int* idx, float **s, float **y, void* ctx);

	struct SQWHParams p; // Model parameters

	float **y;    // Solution vector at grid points y[1..10][1..M]
	float **s;    // Derivative array
	float ***c;   // Auxiliary array for temporary storage

	// Auxiliary arrays
	float *pot;    // Potential array
	float *scalv;  // Scaling array
	int   *indexv; // Index array
	float **z[4];  // Zero field solutions

	// Current parameters
	float    pscale; // Potential scaling factor
	unsigned n;      // Landau level number
	float    H;      // Magnetic field
	float    sH;     // sqrt(H)
	float   **D;     // Derivative matrix

	// Derived parameters
	float mh, ml;
	float A, B;
	float n0, n1, n2, n3;
	float M, N;
	float P, Q;
};
