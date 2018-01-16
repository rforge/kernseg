/*
 *  
 *  Candidate.cpp
 *
 *  Created by Guillem Rigaill on 29/05/13. 
 *  Copyright 2013 INRA, UEVE. All rights reserved.
 *
 */

#include "Candidate.h"
// first line of cumsumXs = 0
double costSegment(double ** cumsumXs, double * cumsumX2, int i, int t, int p){
	double tmp = 0.0;
	double sum;
	for(int j=0; j < p; j++)
		{	
			sum = cumsumXs[i+1][j] - cumsumXs[t+1][j];
			//std::cout << "Sum " << sum << std::endl;
			tmp = tmp + sum*sum;
		}
	//std::cout << "SumX2 " << cumsumX2[i+1] - cumsumX2[t+1] << std::endl;
	return(  cumsumX2[i+1] - cumsumX2[t+1] - tmp / (i -t) );
	

}

