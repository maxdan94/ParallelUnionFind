/*
Maximilien Danisch

Info:
Feel free to use these lines as you wish. 
Efficient implementation of Kruskal's algorithm using a UnionFind datastruture.
- https://en.wikipedia.org/wiki/Kruskal's_algorithm
- https://en.wikipedia.org/wiki/Disjoint-set_data_structure
- https://papers-gamma.link/paper/193

Should scale to at least one billion edges on a commodity machine.

To compile:
"gcc UnionFind.c -O9 -o UnionFind".

To execute:
./UnionFind edgelist.txt
edgelist.txt should contain one edge on each line "u v" u and v are node id (unsigned long long int)
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

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
	unsigned char *r;//ranks
} unionfind;

unionfind* allocuf(unsigned long long n){
	unsigned long long i;
	unionfind* uf=malloc(sizeof(unionfind));
	uf->n=n;
	uf->p=malloc(n*sizeof(unsigned long long));
	uf->r=malloc(n*sizeof(unsigned char));
	for (i=0;i<n;i++){
		uf->p[i]=i;
		uf->r[i]=0;
	}
	return uf;
}

//Find the cluster of element x
unsigned long long Find(unsigned long long x, unionfind *uf){
	if (uf->p[x]!=x){
		uf->p[x]=Find(uf->p[x],uf);
	}
	return uf->p[x];
}

//Merge the clusters xr and yr
void Union(unsigned long long xr, unsigned long long yr, unionfind *uf){
	if (uf->r[xr] < uf->r[yr]){
     		uf->p[xr] = yr;
	}
	else if (uf->r[xr] > uf->r[yr]) {
		uf->p[yr] = xr;
	}
	else {
		uf->p[yr] = xr;
		uf->r[xr] = uf->r[xr]+1;
	}
}

unsigned long long kruskal(edgelist* el){
	unsigned long long i,u,v,p,q,e=0;
	unionfind *uf=allocuf(el->n);

	for (i=0;i<el->e;i++){
		u=el->edges[i].s;
		v=el->edges[i].t;
		p=Find(u,uf);
		q=Find(v,uf);
		if (p!=q){
			Union(p,q,uf);
			e++;
		}
	}
	return e;
}

int main(int argc,char** argv){
	edgelist *el;
	unsigned long long e;

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edge list from file %s\n",argv[1]);
	el=readedgelist(argv[1]);

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
