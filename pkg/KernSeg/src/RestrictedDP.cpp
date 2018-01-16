/*
 *  
 *  RestrictedDP.cpp
 *
 *  Created by Guillem Rigaill on 29/05/13. 
 *  Copyright 2013 INRA, UEVE. All rights reserved.
 *
 */

#include "RestrictedDP.h"


double cost_RestrictedDP(double ** cumsumXs, double * cumsumX2, double *cumweights, int i, int t, int p){
	double tmp = 0.0;
	double sum;
	for(int j=0; j < p; j++)
		{	
			sum = cumsumXs[i+1][j] - cumsumXs[t+1][j];
			//std::cout << "Sum " << sum << std::endl;
			tmp = tmp + sum*sum;
		}
	double weight = cumweights[i+1] - cumweights[t+1];
	//std::cout << "SumX2 " << cumsumX2[i+1] - cumsumX2[t+1] << std::endl;
	return(  cumsumX2[i+1] - cumsumX2[t+1] - tmp / (weight) );
}
	


void Call_RestrictedDP (double * cumsumXs_, double * cumsumX2, double *cumweights, int * n_, int * P_, int * posChanges_, 
		int *nChanges_, int *Kmax_, int * Mkt, double * Ckt){


	int n = *n_;
	int nChanges = *nChanges_;
	int P = *P_;
	int Kmax = *Kmax_;
	int Kmax_tmp;
	double ** cumsumXs = new double*[n+1];
	for(int i=0; i< n+1; i++)
		cumsumXs[i] = new double[P];
	for(int j=0; j < P; j++)
			cumsumXs[0][j] = 0.0;
	for(int i =0; i< n; i++)
		for(int j=0; j < P; j++)
			cumsumXs[i+1][j] = cumsumXs_[i+ n*j];
	
	
	
	for(int t = 0; t < nChanges; t++){
		double test = cost_RestrictedDP(cumsumXs, cumsumX2, cumweights, posChanges_[t], -1, P);
		//std::cout << "posChanges : " << t << ", " <<  posChanges_[t] << " : " << test << std::endl;
		
		Ckt[t] = cost_RestrictedDP(cumsumXs, cumsumX2, cumweights, posChanges_[t], -1, P);
		Mkt[t] = -1;
	}
	
	int tmp, end;
	double * tmpCost = new double[nChanges];
	double newCost;
	
	
	for(int t = 1; t < nChanges;t++){
		// get cost last segment
		//std::cout << "t : " << t << std::endl;
		for(int tau = 0; tau < t; tau++)
			tmpCost[tau] = cost_RestrictedDP(cumsumXs, cumsumX2, cumweights, posChanges_[t], posChanges_[tau], P);
		
		
		if(t < Kmax) {
			Kmax_tmp = t+1;
			} else {
			Kmax_tmp = Kmax;
			}
		//std::cout << "Kmax : " << Kmax  << ", Kmax_tmp : " << Kmax_tmp << std::endl;
		for(int k=2; k <= Kmax_tmp; k++){
			//std::cout << "t : " << t  << ", k : " << k << "/" << Kmax_tmp << std::endl;
			tmp= nChanges*(k-2);
			end= nChanges*(k-1) + t;  

			Ckt[end] = Ckt[tmp+k-2] + tmpCost[k-2]+1;		
			for(int tau = k-2; tau < t; tau++){
				//std::cout << "pos : " << end << ", Curr : " << Ckt[end] << ", prop : " << Ckt[tmp+tau] << " add : "<< tmpCost[tau] << std::endl;
				newCost = Ckt[tmp+tau] + tmpCost[tau];	
				if(newCost < Ckt[end]){
					Ckt[end] = newCost;
					Mkt[end] = tau;
				}
			}
		}
	}
	
	delete[] tmpCost;
	
	for(int i =0; i< n+1; i++)
		delete[] cumsumXs[i];
	delete[] cumsumXs;

}
	
