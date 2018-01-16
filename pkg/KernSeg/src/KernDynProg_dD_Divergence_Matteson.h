#ifndef KernDynProg_dD_DivMatt_h
#define KernDynProg_dD_DivMatt_h

#include <math.h>
#include <stdio.h> 
#include <R.h>
#include <Rinternals.h>

static const double NUMLIB_POSINF=(+1.0/0.0);

#define SQUARE(x) (x)*(x)

void KernDynProg_dD_DivMatt_Raccess(double * vector_mat_i, int *n_, int *d_, int * K_, int *M_ij, double *C_ij, double *delta, double* D_i, int * min_size_, int *option_);


// function for KernSeg like Mattesson divergence (division by n per segment)
void KernDynProg_dD_DivMatt_op0(double ** x_i, int *n_, int *d_, int * K_, int *M_ij, double *C_ij, double *delta, double* D_i, int * min_size_);
// function for Mattesson like Kernseg (division by (n-1) more or less per segment)
// not used anymore equivalent to op0
// void KernDynProg_dD_DivMatt_op1(double ** x_i, int *n_, int *d_, int * K_, int *M_ij, double *C_ij, double *delta, double* D_i, int * min_size_);

#endif // 
