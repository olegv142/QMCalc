#pragma once

#include "errs.h"

// Parameters definitions for single quantum well holes model (sqwh.cpp)
// Units are selected to be dimensionless, see also comments in sqwh.par
struct SQWHParams {
	unsigned M;          // The number of grid points (in z direction)
	unsigned NL;         // Number of cyclotron levels computed
	unsigned subband;    // The quantized subband computed (zero based)
	bool     honly;      // Heavy holes only
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
	void set_params(unsigned n, float field);
	void init_boundary_condition();

	void solve_once(unsigned spin, float precision);
	void solve_zero_field();
	void solve_level(unsigned l);

	unsigned count_zeros(unsigned spin, float precision) const;
	void save_wavefunction(float **f, const std::string& filename) const;
	void save_wavefunctions(unsigned n, unsigned b, const std::string& bname) const;
	void save_levels(const std::string& filename) const;

	static void eq_cb(int k, int* idx, float **s, float **y, void* ctx);
	void eq(int k, int* idx, float **s, float **y) const;
	void get_derivatives(float *a, float **D, int shift, float mult) const;
	void get_equation(float **m, float E) const;
	float barrier_integral(float const fb[8], float E) const;
	void solution_normalize(float **f, float e);

	bool skip_light_hole(unsigned spin) const;

	struct SQWHParams p; // Model parameters

	float **y;      // Solution vector at grid points y[1..10][1..M]
	float **s;      // Derivative array
	float ***c;     // Auxiliary array for temporary storage

	// Auxiliary arrays
	float *scalv;   // Scaling array
	int   *indexv;  // Index array

	// Solutions, the [4] index is spin
	float  **sol_z[4];  // Zero field solutions
	float ****sol_f[4]; // Wavefunctions for different levels and field values
	float **sol_e[4];   // Energy eigenvalues [cyclotron_level][field]

	// Current parameters
	float    H;      // Magnetic field
	float    sH;     // sqrt(H)
	float   **B0;    // Left boundary matrix
	float   **B1;    // Right boundary matrix

	// Derived parameters
	float mh, ml;
	float E0, E1, E2, E3;
	float A, B, M, N, P, Q;
};
