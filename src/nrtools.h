#pragma once

#include "errs.h"

ErrorCode( SOLVDE_IT );

typedef void (*DiffEqCb)( int, int*, float**, float**, void* );

int solvde( DiffEqCb difeq, int itmax, float conv, float slowc, float scalv[], int indexv[],
        int ne, int nb, int m, float **y, float ***c, float **s, FILE *debug, void* ctx );

float pythag(float a, float b);

void gaussj( float **a, int n, float **b, int m );
void gaussj( double **a, int n, double **b, int m );
void ludcmp( float **a, int n, int *indx, int *d );
void lubksb( float **a, int n, int *indx, float *b );
void svdcmp( float **a, int m, int n, float w[], float **v );
void svbksb( float **u, float w[], float **v, int m, int n, float b[], float x[] );

void tred2(float **a, int n, float d[], float e[], bool vect = 1 );
void tqli(float d[], float e[], int n, float **z, bool vect = 1 );

int invset( float **m, int n, int *ind );

