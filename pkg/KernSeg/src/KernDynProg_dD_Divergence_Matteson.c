

#include "KernDynProg_dD_Divergence_Matteson.h"

// Define the appropriate Kernal function


double kernel_dD_DivMatt(double * x, double * y, int d, double * delta){
	double sum1 = 0.0;
   	double sum2 = 0.0;
	double sum3 = 0.0;
	for(int t = 0; t < d; t++){
		sum1 = sum1 + ( x[t] - y[t] ) * ( x[t] - y[t] );
		sum2 = sum2 + x[t] * x[t];
		sum3 = sum3 + y[t] * y[t];
	}
	// alpha = delta[0];
	
	return(pow(sqrt(sum2), delta[0]) + pow(sqrt(sum3), delta[0]) -pow(sqrt(sum1), delta[0]));
}


#define KERN(x, y, delta)  (kernel_dD_DivMatt( (x),  (y), *d_, delta ) )

void KernDynProg_dD_DivMatt_Raccess(double * vector_mat_i, int *n_, int *d_, int * K_, int *M_ij, double *C_ij, double *delta, double* D_i, int * min_size_, int *option_){
	int n = *n_;
	int d = *d_;
	// allocation 
	double** x_i;
	x_i = (double**) malloc( n * sizeof(double*));
	for(int t = 0; t < n ; t++){
		x_i[t] = (double*) malloc( d * sizeof(double));
	}
	
	// fill the vector
	int i_vector_mat = 0;
	for(int t = 0; t < n ; t++){
		for(int i = 0; i < d ; i++){
			//printf("%f, ", vector_mat_i[i_vector_mat]);
			x_i[t][i] = vector_mat_i[i_vector_mat];
			i_vector_mat++;
		}
		//printf("\n");

	}
	if(option_[0] == 0){
	  KernDynProg_dD_DivMatt_op0(x_i, n_, d_, K_, M_ij, C_ij, delta, D_i, min_size_);
	}
	
	// free
	for(int t = 0; t < n; t++){
		free(x_i[t]);
	}
	free(x_i);

}


void KernDynProg_dD_DivMatt_op0(double ** x_i, int *n_, int *d_, int * K_, int *M_ij, double *C_ij, double *delta, double* D_i, int * min_size_){
	int min_size = min_size_[0];
        //#include "Core_1D"

/////////////////////////////////////////////////////////////////////////////////////////////
/// Core of KernSeg BEGIN (copied in all KernDyn....c files)
/////////////////////////////////////////////////////////////////////////////////////////////
	int n = *n_;
	int K = *K_;
	int i, t, k, l, tau, end, tmp;
	double  newCost;
	// ---------------------------------------------- 
	// Initialization
	// ---------------------------------------------- 
	// Initialisation compute K(x_i, x_i) 
	// Depend on the kernel
	//D_i = (double *) malloc( n * sizeof(double));
  	for(t =0; t < n; t++) D_i[t] = KERN(x_i[t], x_i[t], delta);
	

	// free(D_i);

	// Initialisation : l_t and others 
	double * l_t = (double *) malloc( n * sizeof(double));
	double * Sv_t = (double *) malloc( n * sizeof(double));
	double * Cv_t =(double *) malloc( n * sizeof(double));
	for(t =0; t < n; t++) l_t[t] = 0;
	for(t =0; t < n; t++) Sv_t[t] = 0;
	for(t =0; t < n; t++) Cv_t[t] = 0;


	// Initialisation C_ij for i = 0
	// take n^2 
	
	double sum_=0;
	for(t =0; t < n; t++){
		sum_ = sum_ + KERN(x_i[t], x_i[t], delta); 
		for(i =0; i < t; i++){
			 sum_ = sum_ + 2* KERN(x_i[t], x_i[i], delta);
		}
		if( t+1 <= min_size) 
		{
			C_ij[t] = INFINITY ;
		} else 
		{
			C_ij[t] = -sum_ / (t+1);
		}
	}
	// ---------------------------------------------- 
	// Main
	// ---------------------------------------------- 
	
	for(t = 1; t < n; t++){
		// l_t update
		for(i = 0; i < t; i++) {
			l_t[i] = l_t[i] +  KERN(x_i[i], x_i[t], delta);
			//printf("t:%d, i:%d, lt:%f --- \n", t, i, l_t[i]);
		}
		
		
		// segment cost calculation
		Sv_t[t] =  D_i[t];
		Cv_t[t] = - D_i[t];
		//printf("t:%d, i:%d, Sv:%f, Cv:%f\n", t, t, Sv_t[i], Cv_t[i]);
		//printf("t:%d - i:%d - Cv:%f \n", t, t, Cv_t[t]);
		l=1;
		for(i = t-1; i>=0; i--){
			l=l+1;
			Sv_t[i] = Sv_t[i+1] + 2 * l_t[i] + D_i[i];
			Cv_t[i] = - Sv_t[i] / l;
			//printf("t:%d, i:%d, Sv:%f, Cv:%f\n", t, i, Sv_t[i], Cv_t[i]);
			//printf("t:%d - i:%d - Cv:%f \n", t, i, Cv_t[i]);

		}
	
		//position to update n*k+t-1
		for(k = 2; k <= K; k++){
			tmp= n*(k-2);
			end= n*(k-1) + t; 
			C_ij[end]= INFINITY;
			M_ij[end] = -1;// meaning no change for now
			//C_ij[end]=C_ij[tmp+(k-1)-1] + Cv_t[(k-1)]+1;
			for(tau= k-1; tau <= t; tau++){ 
				if(t-tau >= min_size){ // condition for minimum size ?
					newCost = C_ij[tmp+tau-1] + Cv_t[tau];
					//printf("AA # k:%d - tau:%d, %f - %f = cost:%f- bC:%f, end: %d\n", k, tau,  
					//Cv_t[tau], C_ij[tmp+tau-1], newCost, C_ij[end], end);

					if(newCost < C_ij[end]){
						C_ij[end] = newCost;
						M_ij[end] = tau;
					}
				}
				
			}
			
		}
	
	}


	// ---------------------------------------------- 
	// Free
	// ---------------------------------------------- 
	free(l_t);
	free(Sv_t);
	free(Cv_t);

/////////////////////////////////////////////////////////////////////////////////////////////
/// Core of KernSeg END
/////////////////////////////////////////////////////////////////////////////////////////////
	}


// for now does not seems to work properly...
//void KernDynProg_dD_DivMatt_op1(double ** x_i, int *n_, int *d_, int * K_, int *M_ij, double *C_ij, double *delta, double* D_i, int * min_size_){
//	int min_size = min_size_[0];
//	#include "Matteson_1D"
//	}

