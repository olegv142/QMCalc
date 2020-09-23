#pragma once

#include "errs.h"

// Parameters definitions for single quantum well holes model (sqwh.cpp)
// Units are selected to be dimensionless, see also comments in sqwh.par
struct SQWHParams {
	unsigned M;          // The number of grid points (in z direction)
	unsigned NL;         // Number of cyclotron levels computed
	unsigned subband;    // The quantized subband computed (zero based)
	unsigned psteps;     // Number of steps for potential turning on 
	float    hb;         // Barriers height in dimensionless units
	float    Bstep;      // Magnetic field step (in dimensionless units)
	unsigned Bsteps;     // Magnetic field steps (so B varies from 0 to B_step * B_steps)
	unsigned Bskip;      // Write one per Bskip B steps to output file
	float    ef;         // Poisson equation factor (phy" + ef * psy^2 = 0)
	float    me;         // Electrons effective mass
	float    Eo;         // Energy scale in output files
	float    Bo;         // Magnetic field scale in output files
	float    g1, g2, g3, K; // Luttinger parameters
	float    e_side;     // 0 - doping at left side, 1 - doping at right side, .5 - equal doping at both sides
	float    e_cyc;      // Electron cyclotron energy slope (for transition energy calculation)
	float    e_spin;     // Electron spin splitting slope (for transition energy calculation)
	float    prec;       // Solving precision
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
	void init_eguess();
	void init_guess(unsigned spin);
	void set_params(unsigned n, float field);
	void init_boundary_condition();

	void solve_electrons();
	void solve_e(float precision);
	void solve_h(unsigned spin, float precision);
	void solve_zero_field();
	void solve_level(unsigned l);

	void save_wavefunction(float **f, const std::string& filename) const;
	void save_wavefunctions(unsigned n, unsigned b, const std::string& bname) const;
	void save_levels(const std::string& filename, bool transitions) const;

	static void eqe_cb(int k, int* idx, float **s, float **y, void* ctx);
	static void eqh_cb(int k, int* idx, float **s, float **y, void* ctx);
	void  eqe(int k, int* idx, float **s, float **y) const;
	void  eqh(int k, int* idx, float **s, float **y) const;
	void  get_derivatives(float *a, float hpot, float **D, int shift, float mult) const;
	void  get_equation(float **m, float E) const;
	float barrier_integral(float const fb[8], float E) const;
	void  solution_normalize(float **f, float e);
	float get_intensity(float **f, unsigned n0) const;

	struct SQWHParams p; // Model parameters

	float **y;      // Solution vector at grid points y[1..10][1..M]
	float **s;      // Derivative array
	float ***c;     // Auxiliary array for temporary storage
	float *pot;     // Electrostatic potential array

	// Auxiliary arrays
	float *scalv;   // Scaling array
	int   *indexv;  // Index array

	// Solutions, the [4] index is spin
	float *sol_el;     // Electrons wavefunction
	float **sol_z[4];  // Zero field solutions
	float ****sol_f[4];// Wavefunctions for different levels and field values
	float **sol_e[4];  // Energy eigenvalues [cyclotron_level][field]

	// Current parameters
	float    pscale; // Potential scaling factor
	float    H;      // Magnetic field
	float    sH;     // sqrt(H)
	float   **B0;    // Left boundary matrix
	float   **B1;    // Right boundary matrix

	// Derived parameters
	float mh, ml;
	float E0, E1, E2, E3;
	float A, B, M, N, P, Q;
};
