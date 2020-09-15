#include "misc.h"
#include "nrutil.h"
#include "nrtools.h"
#include "math.h"

DefineErrorMess( EMOD_IT, "Too many iterations in emod()" );

/* Count zero crossings given the vector f[1..r] */
unsigned count_zeros(float const *f, int l, int r, float delta)
{
	unsigned cnt = 0;
	int i, sign = 0;
	for( i = l ; i <= r ; i++ ) {
		float val = f[i];
		if ( val > delta ) {
			if ( sign < 0 ) ++cnt;
			sign = 1;
		} else if ( val < -delta ) {
			if ( sign > 0 ) ++cnt;
			sign = -1;
		}
	}
	return cnt;
}

/* Find matrix m[1..n][1..n] with the set of eigenvalues that may differ from eigenvalues of a[1..n][1..n]
 * by the sign of the real part such that all m eigenvalues have positive real part.
 */
void emod(float const * const *a, int n, float conv, int itmax, float **m)
{
	float ** mi = matrix(1, n, 1, n);
	float ** a2 = matrix(1, n, 1, n);
	float sum;
	int i, j, k, it;

	clear_matrix(m, 1, n, 1, n);
	for (i = 1; i <= n; ++i)
	{
		sum = 0;
		for (k = 1; k <= n; ++k)
			sum += fabs(a[i][k] * a[k][i]);
		// Build initial guess
		m[i][i] = sqrt(sum);
		for (j = 1; j <= n; ++j)
		{
			sum = 0;
			for (k = 1; k <= n; ++k)
				sum += a[i][k] * a[k][j];
			// Build a^2 matrix
			a2[i][j] = sum;
		}
	}
	/* Solve a^2 = m^2 equation using Babylonian iteration */
	for (it = 0; it < itmax; ++it)
	{
		float x, err, max_err = 0;
		copy_matrix(m, mi, 1, n, 1, n, 1, 1);
		gaussj(mi, n, NULL, 0);
		for (i = 1; i <= n; ++i)
		{
			for (j = 1; j <= n; ++j) {
				sum = 0;
				for (k = 1; k <= n; ++k)
					sum += a2[i][k] * mi[k][j] + mi[i][k] * a2[k][j];
				x = m[i][j] / 2 + sum / 4;
				err = fabs(x - m[i][j]);
				if (err > max_err)
					max_err = err;
				m[i][j] = x;
			}
		}
		if (max_err < conv)
			break;
	}

	free_matrix(mi, 1, n, 1, n);
	free_matrix(a2, 1, n, 1, n);
	if (it >= itmax)
		Signal(EMOD_IT);
}

