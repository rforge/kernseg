


#include "sumBlocks.h"
double kernelBlock(double x, double y, double delta){
	double sum = (x - y) * (x - y);
	return(exp(-sum / delta));
}


void sumsBlock_Gaussian(double *x_i , int *n_, double *Sums_, double *delta_) {

double delta = delta_[0];
int n = n_[0];
int current_index = 0;

for(int start =0; start < (n-1); start++){
  for(int end = start+1; end < n; end++){

    double tmp = 0;
    for(int i = start; i < end; i++){
      for(int j = start; j <end; j++){
        tmp = tmp + kernelBlock(x_i[i], x_i[j], delta);
      }
    }
    Sums_[current_index] = tmp;
    current_index = current_index+1;

  }
}



}
