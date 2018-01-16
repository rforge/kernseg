/*
 *  Node.h
 * 
 *
 *  Created by Guillem Rigaill on 16/05/11. Modified from Michel Koskas
 *  Copyright 2013 INRA, UEVE. All rights reserved.
 *
 */


#ifndef _Node_h_
#define _Node_h_


class Node
{
public:
	int Index;
	double Value;
	int LowIndex;
	int HighIndex;
	Node(int I, double V, int LI, int HI);
	Node();
	bool operator<(Node &Other);
	bool operator<=(Node &Other);
	Node operator=(Node &Other);
};

void Swap(Node &a, Node &b);


#endif


