/*
 *  
 *  RestrictedDP.h
 *
 *  Created by Guillem Rigaill on 29/05/13. 
 *  Copyright 2013 INRA, UEVE. All rights reserved.
 *
 */

#ifndef _RestrictedDP_
#define _RestrictedDP_

//#include <iostream>

extern "C" {
double cost_RestrictedDP(double ** cumsumXs, double * cumsumX2, double *cumweights, int * posChanges, int i, int t, int p);
void Call_RestrictedDP (double * cumsumXs_, double * cumsumX2, double *cumweights, int * n_, int * P_, int * posChanges_, 
		int *nChanges_, int *Kmax_, int * Mkt, double * Ckt);
}
#endif
