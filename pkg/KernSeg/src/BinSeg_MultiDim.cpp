/*
 *  Binseg.cpp
 *  
 *
 *  Created by Guillem Rigaill on 16/05/11. Modified from Michel Koskas
 *  Copyright 2013 INRA, UEVE. All rights reserved.
 *
 */

#include "BinSeg_MultiDim.h"

BinSeg_MultiDim::BinSeg_MultiDim(double **data_, double * weights_, int n_, int P_, int nbAllocSize_)
{
	data = data_;
	weights = weights_;
	n = n_;
	P= P_;
	//Node *TheN = NULL;
	//MyHeap(nbAllocSize_);
}

/*
double BinSeg_MultiDim::Cost()
{
	int FirstIndex, LastIndex = 0;
	double Res = 0;
	for (std::list<int>::iterator I = Ruptures.begin(); I != Ruptures.end(); I++)
	{
		FirstIndex = LastIndex;
		LastIndex = *I;
		double *Means = new double[P];
		for(int i_p = 0; i_p < P; i_p++)
				Means[i_p] =0;

		for (int i = FirstIndex; i < LastIndex; i++)
			for(int i_p = 0; i_p < P; i_p++)
				Means[i_p] += data[i][i_p];
		
		for(int i_p = 0; i_p < P; i_p++)
			Means[i_p] = Means[i_p]/(LastIndex - FirstIndex);

		for (int i = FirstIndex; i < LastIndex; i++)
			for(int i_p = 0; i_p < P; i_p++)
				Res += (data[i][i_p] - Means[i_p]) * (data[i][i_p] - Means[i_p]);
		delete[] Means;
	}
	return Res;
}*/

Node  BinSeg_MultiDim::Best(int LowIndex, int HighIndex)
{
	if (HighIndex - LowIndex <= 1)
	{
	}
	double SquaresSum = 0;
	double *LeftSum = new double[P];
	for(int i_p = 0; i_p < P; i_p++)
				LeftSum[i_p] =0;
	double *RightSum = new double[P];
	for(int i_p = 0; i_p < P; i_p++)
				RightSum[i_p] =0;

	// Get Sum of Squares for every columns
	for (int i = LowIndex; i < HighIndex; i++)
		for(int i_p = 0; i_p < P; i_p++)
			SquaresSum += weights[i] * (data[i][i_p] * data[i][i_p]); // change W
	//std::cout << "SqSum" << SquaresSum << std::endl;
	double totalWeights=0.0;
	for (int i = LowIndex; i < HighIndex; i++) totalWeights=totalWeights+weights[i];
	
	double leftWeights = 0;
	double rightWeights = totalWeights;
	// Get Sum for each columns
	for (int i = LowIndex; i < HighIndex; i++)
		for(int i_p = 0; i_p < P; i_p++)
				RightSum[i_p] = RightSum[i_p] + weights[i] *data[i][i_p]; // change W
	//for(int i_p = 0; i_p < P; i_p++)
			//	std::cout << i_p << ": " << RightSum[i_p] << std::endl;

	// Get Whole Segment Cost
	double WholeSegmentCost = 0.0; //SquaresSum;
	
	for(int i_p = 0; i_p < P; i_p++)
		WholeSegmentCost =  WholeSegmentCost - RightSum[i_p] * RightSum[i_p]; /// (HighIndex - LowIndex);
	//std::cout << "WholeCost'" << WholeSegmentCost << "SqSum" << SquaresSum << std::endl;
	
	WholeSegmentCost = SquaresSum + WholeSegmentCost/totalWeights; /// (HighIndex - LowIndex) ;
	
	//std::cout << "WholeCost" << WholeSegmentCost << "SqSum" << SquaresSum << std::endl;

	// Initial segment cost for the split higher the without split cost ...
	//double InitialSegmentCost = SquaresSum - WholeSegmentCost;
	//for(int i_p = 0; i_p < P; i_p++)
	//	InitialSegmentCost = InitialSegmentCost - RightSum[i_p] * RightSum[i_p] / (HighIndex - LowIndex -1);
	double InitialSegmentCost = WholeSegmentCost + 1.0;
	// First possible node = doing nothing Cost > Whole Segment 
	
	Node Result(LowIndex, InitialSegmentCost, LowIndex, HighIndex);
	
	double newSegmentCost = 0;
	double tmpL, tmpR;
	
	for (int i = LowIndex + 1; i < HighIndex; i++) // consider all possible split in two
	{
		leftWeights = leftWeights+weights[i-1];
		rightWeights = rightWeights - weights[i-1];
		for(int i_p = 0; i_p < P; i_p++)
				LeftSum[i_p] = LeftSum[i_p] + weights[i-1]*data[i-1][i_p]; 
		for(int i_p = 0; i_p < P; i_p++)
				RightSum[i_p] = RightSum[i_p] - weights[i-1]*data[i-1][i_p]; 
		newSegmentCost = SquaresSum - WholeSegmentCost;
		tmpL = 0;
		tmpR = 0;
		for(int i_p = 0; i_p < P; i_p++) {
			tmpL = tmpL - LeftSum[i_p] * LeftSum[i_p] ;  //  /(i - LowIndex) 
			tmpR = tmpR - RightSum[i_p] * RightSum[i_p]; // /(HighIndex - i);
			}
			
		newSegmentCost = newSegmentCost + tmpL/ (leftWeights) + tmpR/(rightWeights);
							

		Node NewResult(i, newSegmentCost, LowIndex, HighIndex);
		//std::cout << "Split " << i << "New : " << newSegmentCost << " Old : "<< Result.Value << std::endl;
		if (NewResult < Result)
			Result = NewResult;
	}
	return Result;
}

void  BinSeg_MultiDim::Initialize(int NbRuptures)
{
	int LowIndex = 0, HighIndex = n;
	Node FirstNode = Best(LowIndex, HighIndex);
//	std::cout << "FirstNode" << FirstNode.Index <<  ", " <<
//		FirstNode.LowIndex << ", " <<
//		FirstNode.HighIndex << std::endl;
	
	//MyHeap(FirstNode, 1, 2*NbRuptures+10);
	MyHeap.AddNode(FirstNode);

	for (int i = 0; i < NbRuptures; i++)
	{
	//	std::cout << "Round " << i << std::endl;
		Node N = MyHeap.MyHeap[0];
	//	std::cout << "Best Node : " << N.LowIndex << ", " << N.Index << ", " << N.HighIndex << std::endl;
	//	std::cout << "Heap Lgth : " << MyHeap.HeapSize << std::endl;

		Ruptures.push_back(N.Index);
		RupturesCost.push_back(N.Value);

		MyHeap.RemoveHead();
	//	std::cout << "Heap Lgth' : " << MyHeap.HeapSize << std::endl;

		if (N.Index - N.LowIndex > 1)
		{
			Node Left = Best(N.LowIndex, N.Index);
			MyHeap.AddNode(Left);
	//		std::cout << "Left->" << Left.Index << std::endl;

		}
		if (N.HighIndex - N.Index > 1)
		{
			Node Right = Best(N.Index, N.HighIndex);
			MyHeap.AddNode(Right);
	//		std::cout << "Right->" << Right.Index << std::endl;

		}
	}
}





