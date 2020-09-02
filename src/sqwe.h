#pragma once

#include "errs.h"

// Parameters definitions for single quantum well electron model (sqwe.cpp)
// Units are selected to be dimensionless, see also comments in sqwe.par
struct SQWEParams {
	unsigned M;        // The number of grid points (in z direction)
	unsigned subbands; // Number of electron's subbands computed (>0)
	unsigned psteps;   // Number of steps for potential turning on
	float    W;        // SQW width in the units of Z
	float    me;       // Electrons effective mass
	float    me0;      // Electrons effective mass in the barrier (to the left of the well)
	float    me1;      // Electrons effective mass in the barrier (to the right of the well)
	float    ef;       // Poisson equation factor (phy" + ef * psy^2 = 0)
	float    e1;       // External electric field
	float    cex;      // Exchange energy factor (times psy ^ 2/3)
	float    eb0;      // Barriers height for electrons (to the left of the well)
	float    eb1;      // Barriers height for electrons (to the right of the well)
	float    Eo;       // Energy scale in output files
};

class TExpression;

void GetSQWEParams(struct SQWEParams& p, TExpression const& e);

ErrorCode( SQWE_NLS );
ErrorCode( SQWE_OUT );
ErrorCode( SQWE_QNCH );

class SQWESolver {
public:
	SQWESolver(struct SQWEParams const& params);
	~SQWESolver();

	void  Solve(float precision);
	float Energy(unsigned subband) const;
	void  SaveResults(const char* filename) const;
	const float* WaveFunction(unsigned subband) const;
	struct SQWEParams const& Params() const { return p; }

protected:
	void  eguess();
	void  eintro();
	void  esolve(float precision);
	void  adjustEnergy(float delta);
	void  setEnergy(unsigned subband, float val);
	float getEnergy(unsigned subband) const;
	void  checkZeroes(unsigned subband) const;

	void eq0(int k, int* idx, float **s, float **y) const;
	void eq1(int k, int* idx, float **s, float **y) const;

	static void eq0cb(int k, int* idx, float **s, float **y, void* ctx);
	static void eq1cb(int k, int* idx, float **s, float **y, void* ctx);

	struct SQWEParams p; // Model parameters

	float pscale; // Electrostatic potential scaling factor

	float **s;    // Derivative array
	float ***c;   // Axillary array for temporary storage
	float **yg;   // Ground subband:
	// yg[1][1..M] - electrostatic potential
	// yg[2][1..M] - wavefunction density (psy^2) integral
	// yg[3][1..M] - wavefunction
	// yg[4][1..M] - wavefunction derivative
	// yg[5][1..M] - energy eigenvalue

	float ***ys;  // Upper subbands:
	// ys[1..esb-1][1][1..M] - wavefunction density (psy^2) integral
	// ys[1..esb-1][2][1..M] - wavefunction
	// ys[1..esb-1][3][1..M] - wavefunction derivative
	// ys[1..esb-1][4][1..M] - energy eigenvalue

	// Axillary arrays
	float *pot;
	float *scalv0;
	float *scalvs;
	int   *indexv;
};
