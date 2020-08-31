#pragma once

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define nrcheck( cond, text ) if(!(cond)) nrerror( text )
#define nrcheck3( cond, text, value ) if(!(cond)) nrerror( text, value )

void nrerror( char error_text[] );
void nrerror( char error_text[], int value );

float   *vector( int l, int h );
double *dvector( int l, int h );
int    *ivector( int l, int h );
char   *cvector( int l, int h );

float   **matrix( int rl, int rh, int cl, int ch );
double **dmatrix( int rl, int rh, int cl, int ch );
int    **imatrix( int rl, int rh, int cl, int ch );
char   **cmatrix( int rl, int rh, int cl, int ch );

float **submatrix( float **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl );

float  ***f3tensor( int rl, int rh, int cl, int ch, int dl, int dh );
float ****f4tensor( int rl, int rh, int cl, int ch, int dl, int dh, int el, int eh );
float *****f5tensor( 
			int rl, int rh,
			int cl, int ch,
			int dl, int dh,
			int el, int eh,
			int fl, int fh );
float ******f6tensor( 
			int rl, int rh,
			int cl, int ch,
			int dl, int dh,
			int el, int eh,
			int fl, int fh,
			int gl, int gh );

void clear_vector ( float  *v, int l, int h );
void clear_dvector( double *v, int l, int h );
void clear_ivector( int    *v, int l, int h );
void clear_cvector( char   *v, int l, int h );

void clear_matrix ( float  **m, int rl, int rh, int cl, int ch );
void clear_dmatrix( double **m, int rl, int rh, int cl, int ch );
void clear_imatrix( int    **m, int rl, int rh, int cl, int ch );
void clear_cmatrix( char   **m, int rl, int rh, int cl, int ch );

void clear_f3tensor( float  ***t, int rl, int rh, int cl, int ch, int dl, int dh );
void clear_f4tensor( float ****t, int rl, int rh, int cl, int ch, int dl, int dh, int el, int eh );

void copy_vector ( float  *src, float  *dst, int l, int h, int to );
void copy_dvector( double *src, double *dst, int l, int h, int to );
void copy_ivector( int    *src, int    *dst, int l, int h, int to );
void copy_cvector( char   *src, char   *dst, int l, int h, int to );

void copy_matrix ( float  **src, float  **dst, int rl, int rh, int cl, int ch, int tor, int toc );
void copy_dmatrix( double **src, double **dst, int rl, int rh, int cl, int ch, int tor, int toc );
void copy_imatrix( int    **src, int    **dst, int rl, int rh, int cl, int ch, int tor, int toc );
void copy_cmatrix( char   **src, char   **dst, int rl, int rh, int cl, int ch, int tor, int toc );

void copy_vector ( float  *src, float *dst, int l, int h );
void copy_matrix ( float **src, float **dst, int rl, int rh, int cl, int ch );
void copy_f3tensor( float ***src, float ***dst, int rl, int rh, int cl, int ch, int dl, int dh );
void copy_f4tensor( float ****src, float ****dst, int rl, int rh, int cl, int ch, int dl, int dh, int el, int eh );
void copy_f5tensor( float *****src, float *****dst, 
			int rl, int rh,
			int cl, int ch,
			int dl, int dh,
			int el, int eh,
			int fl, int fh );
void copy_f6tensor( float ******src, float ******dst, 
			int rl, int rh,
			int cl, int ch,
			int dl, int dh,
			int el, int eh,
			int fl, int fh,
			int gl, int gh );

void free_vector ( float  *v, int l, int h );
void free_dvector( double *v, int l, int h );
void free_ivector( int    *v, int l, int h );
void free_cvector( char   *v, int l, int h );

void free_matrix ( float  **m, int rl, int rh, int cl, int ch );
void free_dmatrix( double **m, int rl, int rh, int cl, int ch );
void free_imatrix( int    **m, int rl, int rh, int cl, int ch );
void free_cmatrix( char   **m, int rl, int rh, int cl, int ch );

void free_submatrix( float **b, int rl, int rh, int cl, int ch );

void free_f3tensor( float  ***t, int rl, int rh, int cl, int ch, int dl, int dh );
void free_f4tensor( float ****t, int rl, int rh, int cl, int ch, int dl, int dh, int el, int eh );
void free_f5tensor( 	float *****t,
			int rl, int,
			int cl, int,
			int dl, int,
			int el, int,
			int fl, int );
void free_f6tensor( 	float ******t,
			int rl, int,
			int cl, int,
			int dl, int,
			int el, int,
			int fl, int,
			int gl, int );
