#include "misc.h"

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
