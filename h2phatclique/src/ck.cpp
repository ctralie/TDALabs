/*
This is based on code written by:
Maximilien Danisch and Qinna Wang
January 2015
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <algorithm>

#include "ck.h"
//#include "bmatrix.h"

namespace ck
{
bheap *construct(cliques::index n_max){
	cliques::index i;
	bheap* heap = new bheap; 

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=new cliques::index[n_max];
	for (i=0;i<n_max;i++) heap->pt[i]=-1;
	heap->kv=new keyvalue[n_max];
	return heap;
}

inline void swap(bheap *heap,cliques::index i, cliques::index j) {
	keyvalue kv_tmp=heap->kv[i];
	cliques::index pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

inline void bubble_up(bheap *heap,cliques::index i) {
	cliques::index j=(i-1)/2;
	while (i>0) {
		if (heap->kv[j].value>heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

inline void bubble_down(bheap *heap) {
	cliques::index i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->kv[j2].value<heap->kv[j1].value) ) ? j2 : j1 ;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

inline void insert(bheap *heap,keyvalue kv){
	heap->pt[kv.key]=(heap->n)++;
	heap->kv[heap->n-1]=kv;
	bubble_up(heap,heap->n-1);
}

inline void update(bheap *heap,cliques::index key){
	cliques::index i=heap->pt[key];
	if (i!=-1){
		((heap->kv[i]).value)--;
		bubble_up(heap,i);
	}
}

inline keyvalue popmin(bheap *heap){
	keyvalue min=heap->kv[0];
	heap->pt[min.key]=-1;
	heap->kv[0]=heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key]=0;
	bubble_down(heap);
	return min;
}


//compute the maximum of three cliques::index
inline cliques::index max3(cliques::index a,cliques::index b,cliques::index c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//for future use in qsort
int cmpfunc (const void * a, const void * b){
	if (*(cliques::index*)a>*(cliques::index*)b){
		return 1;
	}
	return -1;
}

//Building the graph structure
void mkgraph(sparse *g){
	cliques::index i;
	g->d0=new cliques::index[g->n](); 

	for (i=0;i<g->e;i++) {
		g->d0[g->edges[i].s]++;
		g->d0[g->edges[i].t]++;
	}
	g->cd0=new cliques::index[g->n+1];
	g->cd0[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd0[i]=g->cd0[i-1]+g->d0[i-1];
		g->d0[i-1]=0;
	}

	g->adj0=new cliques::index[2*g->e];

	for (i=0;i<g->e;i++) {
		g->adj0[ g->cd0[g->edges[i].s] + g->d0[ g->edges[i].s ]++ ]=g->edges[i].t;
		g->adj0[ g->cd0[g->edges[i].t] + g->d0[ g->edges[i].t ]++ ]=g->edges[i].s;
	}

	#pragma omp parallel for schedule(dynamic, 1) private(i)
	for (i=0;i<g->n;i++) {
		qsort(&g->adj0[g->cd0[i]],g->d0[i],sizeof(cliques::index),cmpfunc);
	}

}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(sparse *g){
	cliques::index i;
	keyvalue kv;
	bheap* heap=construct(g->n);
	for (i=0;i<g->n;i++){
		kv.key=i;
		kv.value=g->d0[i];
		insert(heap,kv);
	}
	return heap;
}

void freeheap(bheap *heap){
	delete[] heap->pt;
	delete[] heap->kv;
	delete heap;
}

//computing degeneracy ordering and core value
void kcore(sparse* g){
	cliques::index i,j;
	keyvalue kv;
	cliques::index c=0;//the core number
	bheap *heap=mkheap(g);
	g->rank=new cliques::index[g->n];

	for (i=0;i<g->n;i++){
		kv=popmin(heap);
		g->rank[kv.key]=i;
		if (kv.value>c){
			c=kv.value;
		}
		for (j=g->cd0[kv.key];j<g->cd0[kv.key+1];j++){
			update(heap,g->adj0[j]);
		}
	}
	freeheap(heap);
	g->core=c;
}

//Add the special feature to the graph structure: a truncated neighborhood contains only nodes with higher rank
void mkspecial(sparse *g){
	cliques::index i,j,k;
	g->d=new cliques::index[g->n]();
	g->cd=new cliques::index[g->n+1];
	g->adj=new cliques::index[g->e];
	g->cd[0]=0;
	for (i=0;i<g->n;i++) {
		g->cd[i+1]=g->cd[i];
		for (j=g->cd0[i];j<g->cd0[i+1];j++){
			k=g->adj0[j];
			if(g->rank[k]>g->rank[i]){  // {https://www.google.com/#q=Vassilevska+cliques
				g->d[i]++;
				g->adj[g->cd[i+1]++]=k;
			}
		}
	}
}

void freesparse(sparse *g){
	delete[] g->d0;
	delete[] g->cd0;
	delete[] g->adj0;
	delete[] g->rank;
	delete[] g->d;
	delete[] g->cd;
	delete[] g->adj;
	delete g;
}

//store the intersection of list1 and list2 in list3 and return the size of list3 (the 3 lists are sorted)
cliques::index merging(cliques::index *list1, cliques::index s1, cliques::index *list2, cliques::index s2,cliques::index *list3){
	cliques::index i=0,j=0,s3=0;
	cliques::index x=list1[0],y=list2[0];
	while (i<s1 && j<s2){
		if(x<y){
			x=list1[++i];
			continue;
		}
		if(y<x){
			y=list2[++j];
			continue;
		}
		list3[s3++]=x;
		x=list1[++i];
		y=list2[++j];
	}
	return s3;
}

void record_clique(cliques::boundary bdry, bmatrix::cliquearray& highercliques, int MinSize)
{
	cliques::clique c(bdry,1);
	highercliques[bdry.size()-MinSize].emplace_back(c);
}


//the recursion to compute all possible intersections
void recurse(cliques::index kmax, cliques::index k, cliques::index* merge, cliques::index* size, sparse* g, cliques::index* nck, cliques::boundary prev, bmatrix::cliquearray& highercliques, int MinSize, bool count){
        cliques::index t=(k-3)*g->core,t2=t+g->core;
        cliques::index i, u, ct;

        if (k==kmax){
                return;
        }

	cliques::boundary bdry (k+1,0); 
        for (i=0; i<(k-1); i++){
                bdry[i]=prev[i];
        }

	for(i=0; i<size[k-3]; i++){
                u=merge[t+i];
                bdry[k-1]=u;
                prev[k-1]=u;

                size[k-2]=merging(&g->adj[g->cd[u]],g->d[u],&merge[t],size[k-3],&merge[t2]);
                nck[k]+=(cliques::index)size[k-2];

		// On first pass, we should just count the number of simplexes in each dim
		// for pre-allocation.
		if (!count) {
	                // Simplex being considered is prev,u,i for i in (merge[t2+ct] for 0<=ct<size[k-2])
                	for(ct=0; ct<size[k-2]; ct++){
                        	bdry[k]=merge[t2+ct];
				// Boundary must be sorted, but also maintained.
				// So if it isn't sorted, store it before sorting.
				if ( !std::is_sorted(bdry.begin(),bdry.end()) )
				{
					cliques::boundary tempbdry = bdry;
					std::sort(tempbdry.begin(),tempbdry.end());
					record_clique(tempbdry,highercliques,MinSize);
				}
				else
					{ record_clique(bdry,highercliques,MinSize); }
	                }
		}

                recurse(kmax, k+1, merge, size, g, nck, prev, highercliques, MinSize,count);
        }
}


void singlepass(sparse *gr, cliques::index kmax, bmatrix::cliquearray& highercliques, int MinSize, bmatrix::indexvector& dimborders,  bool count){

	cliques::index e,u,v,ct;
	cliques::index *merge,*size;
	cliques::boundary bdry;
	cliques::index *nck=new cliques::index[kmax]();
	cliques::index *nck_p;
	cliques::boundary prev;

	nck[0]=gr->n;
        nck[1]=gr->e;

	if (kmax > 2) {
		merge=new cliques::index[(kmax-2)*gr->core];
                size=new cliques::index[kmax-2];
		nck_p=new cliques::index[kmax]();
		cliques::boundary prev(kmax);

		for(e=0; e<gr->e; e++){
			u=gr->edges[e].s;
                        v=gr->edges[e].t;
			bdry.clear();
			bdry.insert(bdry.end(),{u,v});
			prev[0]=u;
			prev[1]=v;

			// find 3-cliques connected to given edge e
			size[0]=merging(&(gr->adj[gr->cd[u]]),gr->d[u],&(gr->adj[gr->cd[v]]),gr->d[v],merge);

			// Only record cliques if we're in the second step -- in the first we're just counting for preallocation
			if (!count) {
				// for each 3-clique found, add its vertex number to bdry and record the clique.
 				for(ct=0;ct<size[0];ct++){
					bdry.push_back(merge[ct]);
					// This isn't ideal.  Can we maintain sort order as we go instead?
					if ( !std::is_sorted(bdry.begin(),bdry.end()) ) 
					{
						cliques::boundary tempbdry = bdry;
						std::sort(tempbdry.begin(),tempbdry.end());
						record_clique(tempbdry,highercliques,MinSize);
					}
					else { record_clique(bdry,highercliques,MinSize); }
					bdry.pop_back();
				}
			}
			// TODO: check if dimborders and nck_p are identical....I think they are....
			nck_p[2]+=size[0];

			recurse(kmax,3,merge,size,gr,nck_p,prev,highercliques, MinSize, count);
		}

		// We only need to set dimension borders on the counting run
		if (count)
		{
			for(ct=2;ct<kmax;ct++)
				{dimborders.push_back(nck_p[ct]+dimborders[ct-1]);}
		}

		delete[] nck_p;
		delete[] merge;
		delete[] size;
	}
}

} // end namespace
