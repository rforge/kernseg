/*
 *  Heap.h
 *  SmartBacktrack
 *
 *  Created by Guillem Rigaill on 16/05/11. Modified from Michel Koskas
 *  Copyright 2013 INRA, UEVE. All rights reserved.
 *
 */

#include <list>
#include <iostream>
#include <fstream>
#include "Node.h"

#ifndef _Heap_
#define _Heap_


using namespace std;

class Heap
{
public:
  Node *MyHeap;
  int HeapSize;
  //Heap(Node *TheN = NULL, int NbNodes = 0, int AllocationSize = 0);
  Heap(int AllocationSize);
  Heap();
  ~Heap();
  void AddNode(Node N);
  void RemoveHead();
private:
  //void Debug();
	int AllocatedSize;
	void ReAllocate();
};

//ostream & operator<<(ostream &s, const Heap &H);



#endif



