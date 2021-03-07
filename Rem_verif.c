/*
Maximilien Danisch

Info:
Feel free to use these lines as you wish. 
Efficient implementation of Kruskal's algorithm using a UnionFind datastruture (implemented using Rem's algorithm).
- https://en.wikipedia.org/wiki/Kruskal's_algorithm
- https://en.wikipedia.org/wiki/Disjoint-set_data_structure
- https://papers-gamma.link/paper/193

Should scale to at least one billion edges on a commodity machine.

To compile:
"gcc Rem_verif.c -o Rem_verif -fopenmp -O3".

To execute:
./Rem_verif nthreads edgelist.txt
edgelist.txt should contain one edge on each line "u v" u and v are node id (unsigned long long int)
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

typedef struct {
	unsigned long long s;
	unsigned long long t;
} edge;

typedef struct {
	unsigned long long n;//number of nodes
	unsigned long long e;//number of edges
	edge *edges;//list of edges
} edgelist;

//Compute the maximum of three unsigned long long integers.
inline unsigned long long int max3(unsigned long long int a,unsigned long long int b,unsigned long long int c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

edgelist* readedgelist(char* input){
	unsigned long long e1=NLINKS;
	edgelist *g=malloc(sizeof(edgelist));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(input,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%llu %llu\n", &(g->edges[g->e].s), &(g->edges[g->e].t))==2){
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		g->e++;
		if (g->e==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;

	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

// unionfind structure :
typedef struct {
	unsigned long long n;//number of objects
	unsigned long long *p;//parents
} unionfind;

unionfind* allocuf(unsigned long long n){
	unsigned long long i;
	unionfind* uf=malloc(sizeof(unionfind));
	uf->n=n;
	uf->p=malloc(n*sizeof(unsigned long long));
	for (i=0;i<n;i++){
		uf->p[i]=i;
	}
	return uf;
}

//Merge the clusters of x and y returns 1 iff x and y belonged to the same cluster, 0 otherwise.
bool Union(unsigned long long x, unsigned long long y, unionfind *uf, edgelist *elr){
	unsigned long long tmp;

	while (uf->p[x] != uf->p[y]){
		if (uf->p[x]<uf->p[y]){
			if (x==uf->p[x]){
				if (x==uf->p[x]){
					uf->p[x]=uf->p[y];
					elr->edges[elr->e].s=x;
					elr->edges[elr->e++].t=y;
					return 0;
				}
			}
			tmp=uf->p[x];
			uf->p[x]=uf->p[y];
			x=tmp;
		}
		if (uf->p[x]>uf->p[y]){
			if (y==uf->p[y]){
				if (y==uf->p[y]){
					uf->p[y]=uf->p[x];
					elr->edges[elr->e].s=x;
					elr->edges[elr->e++].t=y;
					return 0;
				}
			}
			tmp=uf->p[y];
			uf->p[y]=uf->p[x];
			y=tmp;
		}
	}
	return 1;
}


//Merge the clusters of x and y returns 1 iff x and y belonged to the same cluster, 0 otherwise.
void Union_verif(unsigned long long x, unsigned long long y, unionfind *uf){
	unsigned long long tmp;
	while (uf->p[x] != uf->p[y]){
		if (uf->p[x]<uf->p[y]){
			if (x==uf->p[x]){
				if (x==uf->p[x]){
					uf->p[x]=uf->p[y];
					return;
				}
			}
			tmp=uf->p[x];
			uf->p[x]=uf->p[y];
			x=tmp;
		}
		if (uf->p[x]>uf->p[y]){
			if (y==uf->p[y]){
				if (y==uf->p[y]){
					uf->p[y]=uf->p[x];
					return;
				}
			}
			tmp=uf->p[y];
			uf->p[y]=uf->p[x];
			y=tmp;
		}
	}
}

edgelist *alloctree(unsigned long long n){
	edgelist *el=malloc(sizeof(edgelist));
	el->edges=malloc(n*sizeof(edge));
	el->n=0;
	el->e=0;
	return el;
}

unsigned long long kruskal(edgelist* el){
	unsigned long long i,u,v,e=0;
	unionfind *uf=allocuf(el->n);
	edgelist* elr;

	time_t t1=time(NULL);
	#pragma omp parallel private(u,v,elr) shared(el,uf,e)
	{
		elr=alloctree(el->n);
		#pragma omp for
		for (i=0;i<el->e;i++){
			u=el->edges[i].s;
			v=el->edges[i].t;
			if (Union(u,v,uf,elr)==0){
				#pragma omp atomic
				e++;
			}
		}
		#pragma omp barrier
		#pragma omp critical
		for (i=0;i<elr->e;i++){
			u=elr->edges[i].s;
			v=elr->edges[i].t;
			Union_verif(u,v,uf);
		}
	}

	time_t t2=time(NULL);

	printf("- Time parallel session = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	return e;
}


int main(int argc,char** argv){
	edgelist *el;
	unsigned long long e;

	omp_set_num_threads(atoi(argv[1]));

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edge list from file %s\n",argv[2]);
	el=readedgelist(argv[2]);


	printf("Number of nodes = %llu\n",el->n);
	printf("Number of edges = %llu\n",el->e);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing minimum spaning tree\n");

	e=kruskal(el);

	printf("Number of edges in resulting spaning tree: %llu\n",e);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}
