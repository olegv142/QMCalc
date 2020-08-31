#include "errs.h"

#include <math.h>

void tred2(float **a, int n, float d[], float e[], bool vect )
/*
Housholder reduction of a real, simmetric matrix a[1..n][1..n]. 
On output, a is replaced by the ortogonal matrix Q effecting the 
transformation. d[1..n] returns the diagonal elements of the 
tridiagonal matrix, and e[1..n] the off-diagonal elements, with 
e[1]=0. Several statements, as noted in comments, can be omitted 
if only eigenvalues are to be found, in which case a contains no 
useful information on output. Otherwise they are to be included. 
*/
{
	int l,k,j,i;
	float scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)	
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;	
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					float *jpt = a[j] + 1;
					float *ipt = a[i] + 1;
					float *end = a[i] + j;
					while( ipt <= end )
						g += *jpt++ * *ipt++;
					jpt = a[j+1] + j;
					ipt = a[i] + j + 1;
					end = a[i] + l; 
					for( ; ipt <= end ; jpt += n ) 
						g += *jpt * *ipt++; 
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					float *jpt = a[j] + 1;
					float *ept = e + 1;
					float *ipt = a[i] + 1;
					float *end = a[i] + j;
					while( ipt <= end )
						*jpt++ -= f * *ept++ + g * *ipt++;
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		if( vect ) {
			l=i-1;
			if (d[i])
				for (j=1;j<=l;j++) {
					g=0.0;
					float *jpt = a[1] + j;
					float *ipt = a[i] + 1;
					float *end = a[i] + l;
					for( ; ipt <= end ; jpt += n )
						g += *ipt++ * *jpt;
					jpt = a[1] + j;
					ipt = a[1] + i;
					end = a[l] + i;
					for( ; ipt <= end; ipt += n, jpt += n )
						*jpt -= g * *ipt;
				}
		}
		d[i]=a[i][i];
		if( vect ) {
			a[i][i]=1.0;
			for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
		}
	}
}
