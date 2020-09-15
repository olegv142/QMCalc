#pragma once

#include "errs.h"

/* Count zero crossings given the vector f[1..r] */
unsigned count_zeros(float const *f, int l, int r, float delta);

/* Find matrix m[1..n][1..n] with the set of eigenvalues that may differ from eigenvalues of a[r..r+n-1][c..c+n-1]
 * by the sign of the real part such that all m eigenvalues have positive real part.
 */
void emod(float const * const *a, int r, int c, int n, float conv, int itmax, float **m);

ErrorCode( EMOD_IT );
