#pragma once

#include "errs.h"

// Parameters definitions for single quantum well holes model (sqwh.cpp)
// Units are selected to be dimensionless, see also comments in sqwh.par
struct SQWHParams {
	unsigned M;          // The number of grid points (in z direction)
	unsigned NL;         // Number of cyclotron levels computed
	unsigned subband;    // The quantized subband computed (zero based)
	float    hb;         // Barriers height in dimensionless units
	float    Bstep;      // Magnetic field step (in dimensionless units)
	unsigned Bsteps;     // Magnetic field steps (so B varies from 0 to B_step * B_steps)
	float    Eo;         // Energy scale in output files
	float    Bo;         // Magnetic field scale in output files
	float    g1, g2, g3, K; // Luttinger parameters
	float    prec;     // Solving precision
};

class TExpression;

void GetSQWHParams(struct SQWHParams& p, TExpression const& e);

ErrorCode( SQWH_OUT );
ErrorCode( SQWH_QNCH );
ErrorCode( SQWH_CONV );

class SQWHSolver {
public:
	SQWHSolver(struct SQWHParams const& params);
	~SQWHSolver();

	void  Solve();
	void  SaveResults(const char* basename) const;

protected:
	void init_guess(unsigned spin);
	void set_params(unsigned lvl, float field);
	void init_derivative_matrix();
	void init_boundary_matrix();

	void solve_once(unsigned spin, float precision);
	void solve_zero_field();
	void solve_level(unsigned l);

	unsigned count_zeros(unsigned spin, float precision) const;
	void save_wavefunction(float **f, const std::string& filename) const;
	void save_levels(const std::string& filename) const;

	void eq(int k, int* idx, float **s, float **y) const;
	static void eq_cb(int k, int* idx, float **s, float **y, void* ctx);

	struct SQWHParams p; // Model parameters

	float **y;      // Solution vector at grid points y[1..10][1..M]
	float **s;      // Derivative array
	float ***c;     // Auxiliary array for temporary storage

	// Auxiliary arrays
	float *scalv;   // Scaling array
	int   *indexv;  // Index array

	float  **sol_z[4];   // Zero field solutions
	float ***sol_h[4];  // High field solutions
	float ***sol_e;     // Energy results [cyclotron_level][field][spin]

	// Current parameters
	float    pscale; // Potential scaling factor
	unsigned n;      // Landau level number
	float    H;      // Magnetic field
	float    sH;     // sqrt(H)
	float   **D;     // Derivative matrix
	float   **B0;    // Left boundary matrix
	float   **B1;    // Right boundary matrix

	// Derived parameters
	float mh, ml;
	float A, B;
	float n0, n1, n2, n3;
	float M, N;
	float P, Q, P2, Q2;
};
