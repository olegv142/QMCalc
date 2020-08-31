#include "nrtools.h"

#include <math.h>

float pythag( float a, float b )
// Computes ( a^2 + b^2 )^(1/2) without destructive underflow or overflow.
{
	float absa = fabs( a );
	float absb = fabs( b );
	if( absa > absb ) {
		float rat = absb / absa;
		return absa * sqrt( 1 + rat * rat );
	} else if( absb ) {
		float rat = absa / absb;
		return absb * sqrt( 1 + rat * rat );
	} else
		return 0;
}


