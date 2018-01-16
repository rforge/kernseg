#ifndef KernDynProg_dD_Gaussian_h
#define KernDynProg_dD_Gaussian_h

#include <math.h>
#include <stdio.h> 
#include <R.h>
#include <Rinternals.h>



#define SQUARE(x) (x)*(x)

void KernDynProg_dD_Gaussian_Raccess(double * vector_mat_i, int *n_, int *d_, int * K_, int *M_ij, double *C_ij, double *delta, double* D_i, int * min_size_);


void KernDynProg_dD_Gaussian(double ** x_i, int *n_, int *d_, int * K_, int *M_ij, double *C_ij, double *delta, double* D_i, int * min_size_);

#endif // 

