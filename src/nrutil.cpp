#include "nrutil.h"
#include "errs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#include "compat.h"

#define NR_END 1

DefineErrorMess( NRC_ERROR, "NRC Error:" );

void nrerror( char error_text[] )
/* Numerical Recipes standart error handler */
{
	Signal( NRC_ERROR, error_text );
}

void nrerror( char error_text[], int value )
/* Numerical Recipes standart error handler */
{
	unsigned const buff_sz = 128;
	char buf[buff_sz+1];
	buf[buff_sz] = 0;
	snprintf( buf, buff_sz, "%u bytes %s", error_text);
	Signal( NRC_ERROR, buf );
}

float   *vector( int l, int h )
/* allocate a float vector with subscript range v[l..h] */
{
	int size;
	float *v = (float*)malloc(size=(h-l+1+NR_END)*sizeof(float));
	nrcheck3( v, "allocation failure in vector()", size );
	return v + NR_END - l;
}

double *dvector( int l, int h )
/* allocate a double vector with subscript range v[l..h] */
{
	int size;
	double *v = (double*)malloc(size=(h-l+1+NR_END)*sizeof(double));
	nrcheck3( v, "allocation failure in dvector()", size );
	return v + NR_END - l;
}

int    *ivector( int l, int h )
/* allocate an int vector with subscript range v[l..h] */
{
	int size;
	int *v = (int*)malloc(size=(h-l+1+NR_END)*sizeof(int));
	nrcheck3( v, "allocation failure in ivector()", size );
	return v + NR_END - l;
}

char *cvector( int l, int h )
/* allocate a char vector with subscript range v[l..h] */
{
	int size;
	char *v = (char*)malloc(size=(h-l+1+NR_END)*sizeof(char));
	nrcheck3( v, "allocation failure in cvector()", size );
	return v + NR_END - l;
}

float   **matrix( int rl, int rh, int cl, int ch )
/* allocate a float matrix with subscruipt range m[rl..rh][cl..ch] */
{
	int size;
	int nrow = rh - rl + 1;
	int ncol = ch - cl + 1;
	float **m = (float**)malloc(size=(nrow+NR_END)*sizeof(float*));
	nrcheck3( m, "allocation failure 1 in matrix()", size );
	m += NR_END;
	m -= rl;
	m[rl] = (float*)malloc(size=(nrow*ncol+NR_END)*sizeof(float));
	nrcheck3( m[rl], "allocation failure 2 in matrix()", size );
	m[rl] += NR_END;
	m[rl] -= cl;
	for( int i = rl + 1 ; i <= rh ; i++ ) m[i] = m[i-1] + ncol;
	return m;
}

double **dmatrix( int rl, int rh, int cl, int ch )
/* allocate a double matrix with subscruipt range m[rl..rh][cl..ch] */
{
	int size;
	int nrow = rh - rl + 1;
	int ncol = ch - cl + 1;
	double **m = (double**)malloc(size=(nrow+NR_END)*sizeof(double*));
	nrcheck3( m, "allocation failure 1 in dmatrix()", size );
	m += NR_END;
	m -= rl;
	m[rl] = (double*)malloc(size=(nrow*ncol+NR_END)*sizeof(double));
	nrcheck3( m[rl], "allocation failure 2 in dmatrix()", size );
	m[rl] += NR_END;
	m[rl] -= cl;
	for( int i = rl + 1 ; i <= rh ; i++ ) m[i] = m[i-1] + ncol;
	return m;
}

int    **imatrix( int rl, int rh, int cl, int ch )
/* allocate an int matrix with subscruipt range m[rl..rh][cl..ch] */
{
	int size;
	int nrow = rh - rl + 1;
	int ncol = ch - cl + 1;
	int **m = (int**)malloc(size=(nrow+NR_END)*sizeof(int*));
	nrcheck3( m, "allocation failure 1 in imatrix()", size );
	m += NR_END;
	m -= rl;
	m[rl] = (int*)malloc(size=(nrow*ncol+NR_END)*sizeof(int));
	nrcheck3( m[rl], "allocation failure 2 in imatrix()", size );
	m[rl] += NR_END;
	m[rl] -= cl;
	for( int i = rl + 1 ; i <= rh ; i++ ) m[i] = m[i-1] + ncol;
	return m;
}

char   **cmatrix( int rl, int rh, int cl, int ch )
/* allocate a char matrix with subscruipt range m[rl..rh][cl..ch] */
{
	int size;
	int nrow = rh - rl + 1;
	int ncol = ch - cl + 1;
	char **m = (char**)malloc(size=(nrow+NR_END)*sizeof(char*));
	nrcheck3( m, "allocation failure 1 in cmatrix()", size );
	m += NR_END;
	m -= rl;
	m[rl] = (char*)malloc(size=(nrow*ncol+NR_END)*sizeof(char));
	nrcheck3( m[rl], "allocation failure 2 in cmatrix()", size );
	m[rl] += NR_END;
	m[rl] -= cl;
	for( int i = rl + 1 ; i <= rh ; i++ ) m[i] = m[i-1] + ncol;
	return m;
}

float **submatrix( float **a, int oldrl, int oldrh, int oldcl, int, int newrl, int newcl )
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	int size;
	int nrow = oldrh - oldrl + 1;
	int ncol = oldcl - newcl;
	float **m = (float**)malloc(size=(nrow+NR_END)*sizeof(float*));
	nrcheck3( m, "allocation failure in submatrix()", size );
	m += NR_END;
	m -= newrl;
	for( int i = oldrl, j = newrl ; i <= oldrh ; i++, j++ ) m[j] = a[i] + ncol;
	return m;
}

float ***f3tensor( int rl, int rh, int cl, int ch, int dl, int dh )
/* allocate a float 3tensor with range t[rl..rh][cl..ch][dl..dh] */
{
	int size;
	int l2 = ch - cl + 1;
	int l3 = dh - dl + 1;

	int s1 = rh - rl + 1;
	int s2 = s1 * l2;
	int s3 = s2 * l3;

	float ***p1 = (float***)malloc(size=(s1+NR_END)*sizeof(float**));
	nrcheck3( p1, "allocation failure 1 in f3tensor()", size );
	p1 += NR_END;

	float **p2 = (float**)malloc(size=(s2+NR_END)*sizeof(float*));
	nrcheck3( p2, "allocation failure 2 in f3tensor()", size );
	p2 += NR_END;

	float *p3 = (float*)malloc(size=(s3+NR_END)*sizeof(float));
	nrcheck3( p3, "allocation failure 3 in f3tensor()", size );
	p3 += NR_END;

	float ***f1 = p1 - rl;
	float **f2 = *p1 = p2 - cl;
		int i;
	for( i = 1 ; i < s1 ; i++ )
		*(++p1) = f2 += l2;
	float *f3 = *p2 = p3 - dl;
	for( i = 1 ; i < s2 ; i++ )
		*(++p2) = f3 += l3;

	return f1;
}

float ****f4tensor( int rl, int rh, int cl, int ch, int dl, int dh, int el, int eh )
{
	int size;
	int l2 = ch - cl + 1;
	int l3 = dh - dl + 1;
	int l4 = eh - el + 1;

	int s1 = rh - rl + 1;
	int s2 = s1 * l2;
	int s3 = s2 * l3;
	int s4 = s3 * l4;

	float ****p1 = (float****)malloc(size=(s1+NR_END)*sizeof(float***));
	nrcheck3( p1, "allocation failure 1 in f4tensor()", size );
	p1 += NR_END;

	float ***p2 = (float***)malloc(size=(s2+NR_END)*sizeof(float**));
	nrcheck3( p2, "allocation failure 2 in f4tensor()", size );
	p2 += NR_END;

	float **p3 = (float**)malloc(size=(s3+NR_END)*sizeof(float*));
	nrcheck3( p3, "allocation failure 3 in f4tensor()", size );
	p3 += NR_END;

	float *p4 = (float*)malloc(size=(s4+NR_END)*sizeof(float));
	nrcheck3( p4, "allocation failure 4 in f4tensor()", size );
	p4 += NR_END;

	float ****f1 = p1 - rl;
	float ***f2 = *p1 = p2 - cl;
	int i;
	for( i = 1 ; i < s1 ; i++ )
		*(++p1) = f2 += l2;
	float **f3 = *p2 = p3 - dl;
	for( i = 1 ; i < s2 ; i++ )
		*(++p2) = f3 += l3;
	float *f4 = *p3 = p4 - el;
	for( i = 1 ; i < s3 ; i++ )
		*(++p3) = f4 += l4;

	return f1;
}

float *****f5tensor( 	int rl, int rh, 
						int cl, int ch, 
						int dl, int dh, 
						int el, int eh,
						int fl, int fh )
{
	int size;
	int l2 = ch - cl + 1;
	int l3 = dh - dl + 1;
	int l4 = eh - el + 1;
	int l5 = fh - fl + 1;

	int s1 = rh - rl + 1;
	int s2 = s1 * l2;
	int s3 = s2 * l3;
	int s4 = s3 * l4;
	int s5 = s4 * l5;

	float *****p1 = (float*****)malloc(size=(s1+NR_END)*sizeof(float****));
	nrcheck3( p1, "allocation failure 1 in f5tensor()", size );
	p1 += NR_END;

	float ****p2 = (float****)malloc(size=(s2+NR_END)*sizeof(float***));
	nrcheck3( p2, "allocation failure 2 in f5tensor()", size );
	p2 += NR_END;

	float ***p3 = (float***)malloc(size=(s3+NR_END)*sizeof(float**));
	nrcheck3( p3, "allocation failure 3 in f5tensor()", size );
	p3 += NR_END;

	float **p4 = (float**)malloc(size=(s4+NR_END)*sizeof(float*));
	nrcheck3( p4, "allocation failure 4 in f5tensor()", size );
	p4 += NR_END;

	float *p5 = (float*)malloc(size=(s5+NR_END)*sizeof(float));
	nrcheck3( p5, "allocation failure 5 in f5tensor()", size );
	p5 += NR_END;

	float *****f1 = p1 - rl;
	float ****f2 = *p1 = p2 - cl;
	int i;
	for( i = 1 ; i < s1 ; i++ )
		*(++p1) = f2 += l2;

	float ***f3 = *p2 = p3 - dl;
	for( i = 1 ; i < s2 ; i++ )
		*(++p2) = f3 += l3;

	float **f4 = *p3 = p4 - el;
	for( i = 1 ; i < s3 ; i++ )
		*(++p3) = f4 += l4;

	float *f5 = *p4 = p5 - fl;
	for( i = 1 ; i < s4 ; i++ )
		*(++p4) = f5 += l5;

	return f1;
}

float ******f6tensor( 	int rl, int rh, 
						int cl, int ch, 
						int dl, int dh, 
						int el, int eh,
						int fl, int fh,
						int gl, int gh )
{
	int size;
	int l2 = ch - cl + 1;
	int l3 = dh - dl + 1;
	int l4 = eh - el + 1;
	int l5 = fh - fl + 1;
	int l6 = gh - gl + 1;

	int s1 = rh - rl + 1;
	int s2 = s1 * l2;
	int s3 = s2 * l3;
	int s4 = s3 * l4;
	int s5 = s4 * l5;
	int s6 = s5 * l6;

	float ******p1 = (float******)malloc(size=(s1+NR_END)*sizeof(float*****));
	nrcheck3( p1, "allocation failure 1 in f6tensor()", size );
	p1 += NR_END;

	float *****p2 = (float*****)malloc(size=(s2+NR_END)*sizeof(float****));
	nrcheck3( p2, "allocation failure 2 in f6tensor()", size );
	p2 += NR_END;

	float ****p3 = (float****)malloc(size=(s3+NR_END)*sizeof(float***));
	nrcheck3( p3, "allocation failure 3 in f6tensor()", size );
	p3 += NR_END;

	float ***p4 = (float***)malloc(size=(s4+NR_END)*sizeof(float**));
	nrcheck3( p4, "allocation failure 4 in f6tensor()", size );
	p4 += NR_END;

	float **p5 = (float**)malloc(size=(s5+NR_END)*sizeof(float*));
	nrcheck3( p5, "allocation failure 5 in f6tensor()", size );
	p5 += NR_END;

	float *p6 = (float*)malloc(size=(s6+NR_END)*sizeof(float));
	nrcheck3( p6, "allocation failure 6 in f6tensor()", size );
	p6 += NR_END;


	float ******f1 = p1 - rl;
	float *****f2 = *p1 = p2 - cl;
	int i;
	for( i = 1 ; i < s1 ; i++ )
		*(++p1) = f2 += l2;

	float ****f3 = *p2 = p3 - dl;
	for( i = 1 ; i < s2 ; i++ )
		*(++p2) = f3 += l3;

	float ***f4 = *p3 = p4 - el;
	for( i = 1 ; i < s3 ; i++ )
		*(++p3) = f4 += l4;

	float **f5 = *p4 = p5 - fl;
	for( i = 1 ; i < s4 ; i++ )
		*(++p4) = f5 += l5;

	float *f6 = *p5 = p6 - gl;
	for( i = 1 ; i < s5 ; i++ )
		*(++p5) = f6 += l6;

	return f1;
}

void clear_vector ( float *v, int l, int h )
{
	assert( v );
	memset( v + l, 0,
		(size_t)((h-l+1)*sizeof(float)) );
}

void clear_dvector( double *v, int l, int h )
{
	assert( v );
	memset( v + l, 0,
		(size_t)((h-l+1)*sizeof(double)) );
}

void clear_ivector( int *v, int l, int h )
{
	assert( v );
	memset( v + l, 0,
		(size_t)((h-l+1)*sizeof(int)) );
}

void clear_cvector( char *v, int l, int h )
{
	assert( v );
	memset( v + l, 0,
		(size_t)((h-l+1)*sizeof(char)) );
}

void clear_matrix ( float **m, int rl, int rh, int cl, int ch )
{
	assert( m );
	int ncol = ch - cl + 1;
	for( ; rl <= rh ; rl++ )
		memset( m[rl] + cl, 0, (size_t)((ncol)*sizeof(float)) );
}

void clear_dmatrix( double **m, int rl, int rh, int cl, int ch )
{
	assert( m );
	int ncol = ch - cl + 1;
	for( ; rl <= rh ; rl++ )
		memset( m[rl] + cl, 0, (size_t)((ncol)*sizeof(double)) );
}

void clear_imatrix( int **m, int rl, int rh, int cl, int ch )
{
	assert( m );
	int ncol = ch - cl + 1;
	for( ; rl <= rh ; rl++ )
		memset( m[rl] + cl, 0, (size_t)((ncol)*sizeof(int)) );
}

void clear_cmatrix ( char **m, int rl, int rh, int cl, int ch )
{
	assert( m );
	int ncol = ch - cl + 1;
	for( ; rl <= rh ; rl++ )
		memset( m[rl] + cl, 0, (size_t)((ncol)*sizeof(char)) );
}

void clear_f3tensor( float ***t, int rl, int rh, int cl, int ch, int dl, int dh )
{
	assert( t );
	int ndep = dh - dl + 1;
	for( ; rl <= rh ; rl++ )
		for( int j = cl ; j <= ch ; j++ )
			memset( t[rl][j] + dl, 0, (size_t)((ndep)*sizeof(float)) );
}

void clear_f4tensor( float ****t, int rl, int rh, int cl, int ch, int dl, int dh, int el, int eh )
{
	assert( t );
	int ne = eh - el + 1;
	for( ; rl <= rh ; rl++ )
		for( int j = cl ; j <= ch ; j++ )
			for( int k = dl ; k <= dh ; k++ )
				memset( t[rl][j][k] + el, 0, (size_t)((ne)*sizeof(float)) );
}

void copy_vector ( float *src, float *dst, int l, int h, int to )
{
	assert( src );
	assert( dst );
	memcpy( dst + to, src + l, (size_t)((h-l+1)*sizeof(float)) );
}

void copy_dvector ( double *src, double *dst, int l, int h, int to )
{
	assert( src );
	assert( dst );
	memcpy( dst + to, src + l, (size_t)((h-l+1)*sizeof(double)) );
}

void copy_ivector ( int *src, int *dst, int l, int h, int to )
{
	assert( src );
	assert( dst );
	memcpy( dst + to, src + l, (size_t)((h-l+1)*sizeof(int)) );
}

void copy_cvector ( char *src, char *dst, int l, int h, int to )
{
	assert( src );
	assert( dst );
	memcpy( dst + to, src + l, (size_t)((h-l+1)*sizeof(char)) );
}

void copy_matrix ( float **src, float **dst,
			int rl, int rh, int cl, int ch, int tor, int toc )
{
	assert( src );
	assert( dst );
	int ncol = ch - cl + 1;
	for( int i = tor ; rl <= rh ; i++, rl++ )
		memcpy( dst[i] + toc, src[rl] + cl, (size_t)((ncol)*sizeof(float)) );
}

void copy_dmatrix ( double **src, double **dst,
			int rl, int rh, int cl, int ch, int tor, int toc )
{
	assert( src );
	assert( dst );
	int ncol = ch - cl + 1;
	for( int i = tor ; rl <= rh ; i++, rl++ )
		memcpy( dst[i] + toc, src[rl] + cl, (size_t)((ncol)*sizeof(double)) );
}

void copy_imatrix ( int **src, int **dst,
			int rl, int rh, int cl, int ch, int tor, int toc )
{
	assert( src );
	assert( dst );
	int ncol = ch - cl + 1;
	for( int i = tor ; rl <= rh ; i++, rl++ )
		memcpy( dst[i] + toc, src[rl] + cl, (size_t)((ncol)*sizeof(int)) );
}

void copy_cmatrix ( char **src, char **dst,
			int rl, int rh, int cl, int ch, int tor, int toc )
{
	assert( src );
	assert( dst );
	int ncol = ch - cl + 1;
	for( int i = tor ; rl <= rh ; i++, rl++ )
		memcpy( dst[i] + toc, src[rl] + cl, (size_t)((ncol)*sizeof(char)) );
}

void copy_vector( float *src, float *dst, int l, int h )
{
	assert( src && dst );
	memcpy( dst + l, src + l, (size_t)((h - l + 1)*sizeof(float)) );
}

void copy_matrix ( float **src, float **dst, int rl, int rh, int cl, int ch )
{
	assert( src && dst );
	memcpy( dst[rl] + cl, src[rl] + cl, 
		(size_t)((rh - rl + 1)*(ch - cl + 1)*sizeof(float)) );
}

void copy_f3tensor( float ***src, float ***dst, int rl, int rh, int cl, int ch, int dl, int dh )
{
	assert( src && dst );
	memcpy( dst[rl][cl] + dl, src[rl][cl] + dl, 
		(size_t)((rh - rl + 1)*(ch - cl + 1)*(dh - dl + 1)*sizeof(float)) );
}

void copy_f4tensor( float ****src, float ****dst, int rl, int rh, int cl, int ch, int dl, int dh, int el, int eh )
{
	assert( src && dst );
	memcpy( dst[rl][cl][dl] + el, src[rl][cl][dl] + el, 
		(size_t)((rh - rl + 1)*(ch - cl + 1)*(dh - dl + 1)*
				 (eh - el + 1)*sizeof(float)) );
}

void copy_f5tensor( float *****src, float *****dst, 
			int rl, int rh,
			int cl, int ch,
			int dl, int dh,
			int el, int eh,
			int fl, int fh )
{
	assert( src && dst );
	memcpy( dst[rl][cl][dl][el] + fl, src[rl][cl][dl][el] + fl, 
		(size_t)((rh - rl + 1)*(ch - cl + 1)*(dh - dl + 1)*
				 (eh - el + 1)*(fh - fl + 1)*sizeof(float)) );
}

void copy_f6tensor( float ******src, float ******dst,
			int rl, int rh,
			int cl, int ch,
			int dl, int dh,
			int el, int eh,
			int fl, int fh,
			int gl, int gh )
{
	assert( src && dst );
	memcpy( dst[rl][cl][dl][el][fl] + gl, src[rl][cl][dl][el][fl] + gl,
		(size_t)((rh - rl + 1)*(ch - cl + 1)*(dh - dl + 1)*
				 (eh - el + 1)*(fh - fl + 1)*(gh - gl + 1)*sizeof(float)) );
}

void free_vector ( float *v, int l, int )
/* free a float vector allocated with vector() */
{
	assert( v );
	free( v + l - NR_END );
}

void free_dvector( double *v, int l, int )
/* free a double vector allocated with dvector() */
{
	assert( v );
	free( v + l - NR_END );
}

void free_ivector( int *v, int l, int )
/* free an int vector allocated with ivector() */
{
	assert( v );
	free( v + l - NR_END );
}

void free_cvector( char *v, int l, int )
/* free a char vector allocated with cvector() */
{
	assert( v );
	free( v + l - NR_END );
}

void free_matrix ( float **m, int rl, int, int cl, int )
/* free a float matrix allocated by matrix() */
{
	assert( m );
	free( m[rl] + cl - NR_END );
	free( m + rl - NR_END );
}

void free_dmatrix( double **m, int rl, int, int cl, int )
/* free a double matrix allocated by dmatrix() */
{
	assert( m );
	free( m[rl] + cl - NR_END );
	free( m + rl - NR_END );
}

void free_imatrix( int **m, int rl, int, int cl, int )
/* free an int matrix allocated by imatrix() */
{
	assert( m );
	free( m[rl] + cl - NR_END );
	free( m + rl - NR_END );
}

void free_cmatrix ( char **m, int rl, int, int cl, int )
/* free a char matrix allocated by cmatrix() */
{
	assert( m );
	free( m[rl] + cl - NR_END );
	free( m + rl - NR_END );
}

void free_submatrix( float **b, int rl, int, int, int )
/* free a submatrix allocated by submatrix() */
{
	assert( b );
	free( b + rl - NR_END );
}

void free_f3tensor( float ***t, int rl, int, int cl, int, int dl, int )
/* free a float f3tensor allocated by f3tensor() */
{
	assert( t );
	free( t[rl][cl] + dl - NR_END );
	free( t[rl] + cl - NR_END );
	free( t + rl - NR_END );
}

void free_f4tensor( float ****t, int rl, int, int cl, int, int dl, int, int el, int )
{
	assert( t );
	free( t[rl][cl][dl] + el - NR_END );
	free( t[rl][cl] + dl - NR_END );
	free( t[rl] + cl - NR_END );
	free( t + rl - NR_END );
}

void free_f5tensor( 	float *****t,
						int rl, int,
						int cl, int,
						int dl, int,
						int el, int,
						int fl, int )
{
	assert( t );
	free( t[rl][cl][dl][el] + fl - NR_END );
	free( t[rl][cl][dl] + el - NR_END );
	free( t[rl][cl] + dl - NR_END );
	free( t[rl] + cl - NR_END );
	free( t + rl - NR_END );
}

void free_f6tensor( 	float ******t,
						int rl, int,
						int cl, int,
						int dl, int,
						int el, int,
						int fl, int,
						int gl, int )
{
	assert( t );
	free( t[rl][cl][dl][el][fl] + gl - NR_END );
	free( t[rl][cl][dl][el] + fl - NR_END );
	free( t[rl][cl][dl] + el - NR_END );
	free( t[rl][cl] + dl - NR_END );
	free( t[rl] + cl - NR_END );
	free( t + rl - NR_END );
}

