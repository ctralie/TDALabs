/*
This is based on code written by:
Maximilien Danisch and Qinna Wang
January 2015
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr
*/
#pragma once

#ifndef _CK
	#define _CK
#endif

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "clique.h"

#include <bmatrix.h>

namespace ck
{

// heap data structure :

typedef struct {
	cliques::index key;
	cliques::index value;
} keyvalue;

typedef struct {
	cliques::index n_max;	// max number of nodes.
	cliques::index n;	// number of nodes.
	cliques::index *pt;	// pointers to nodes.
	keyvalue *kv; 		// nodes.
} bheap;

extern bheap *construct(cliques::index n_max);

extern void swap(bheap *, cliques::index, cliques::index);

extern void bubble_up(bheap *, cliques::index);

extern void bubble_down(bheap *);

extern void insert(bheap *, keyvalue);

extern void update(bheap *, cliques::index); 

extern keyvalue popmin(bheap *);

// graph datastructure:

typedef struct {
	cliques::index s;
	cliques::index t;
	cliques::weight_t weight;
} edge;

typedef struct {
	cliques::index n; //number of nodes
	cliques::index e; //number of edges
	std::vector<edge> edges; //list of edges

	cliques::index *d0; //degrees
	cliques::index *cd0; //cumulative degree: (start with 0) length=dim+1
	cliques::index *adj0; //list of neighbors

	cliques::index *rank; //degeneracy rankings of nodes

	cliques::index *d; //truncated degrees
	cliques::index *cd; //cumulative degree: (start with 0) length=dim+1
	cliques::index *adj; //list of neighbors with higher rank

	cliques::index core; //core number of the graph
} sparse;


//compute the maximum of three cliques::index
extern cliques::index max3(cliques::index,cliques::index,cliques::index);

//reading the edgelist from file
extern sparse* readedgelist(char* edgelist);

//for future use in qsort
extern int cmpfunc (const void * a, const void * b);

//Building the graph structure
extern void mkgraph(sparse *g);

//Building the heap structure with (key,value)=(node,degree) for each node
extern bheap* mkheap(sparse *g);

extern void freeheap(bheap *heap);

//computing degeneracy ordering and core value
extern void kcore(sparse* g);

//Add the special feature to the graph structure: a truncated neighborhood contains only nodes with higher rank
extern void mkspecial(sparse *g);

extern void freesparse(sparse *g);

//store the intersection of list1 and list2 in list3 and return the size of list3 (the 3 lists are sorted)
cliques::index merging(cliques::index *, cliques::index, cliques::index *, cliques::index, cliques::index *);

// record a single clique in the boundary matrix
extern void record_clique(cliques::boundary bdry, bmatrix::cliquearray& highercliques);

//the recursion to compute all possible intersections
extern void recurse(cliques::index kmax, cliques::index k, cliques::index* merge, cliques::index* size, sparse* g, cliques::index* nck, cliques::boundary prev, bmatrix::cliquearray& highercliques, int MinSize, bool count);

//one pass over all k-cliques
extern void singlepass(sparse *g,cliques::index kmax, bmatrix::cliquearray& highercliques, int MinSize, std::vector<cliques::index>& dimborders,  bool count);

sparse* new_read_extdimacs_file(std::istream& dimacs, std::string& dimacs_err, bmatrix::cliquearray& cliquevector, std::vector<cliques::index>&);

} // end namespace
