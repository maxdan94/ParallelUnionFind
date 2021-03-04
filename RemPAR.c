/*
Maximilien Danisch

Info:
Feel free to use these lines as you wish. 
Efficient implementation of Kruskal's algorithm using a UnionFind datastruture (implemented using Rem'a algorithms).
- https://en.wikipedia.org/wiki/Kruskal's_algorithm
- https://en.wikipedia.org/wiki/Disjoint-set_data_structure
- https://papers-gamma.link/paper/193

Should scale to at least one billion edges on a commodity machine.

To compile:
"gcc kruskal.c -O9 -o kruskal".

To execute:
./kruskal edgelist.txt res.txt
edgelist.txt should contain one edge on each line "u v w" u and v are node id (unsigned int) and w is the edge weight (double).
res.txt will contained the list of the edges of the resulting tree
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

typedef struct {
	unsigned s;
	unsigned t;
	//double w;
} edge;

typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges
} edgelist;

//Compute the maximum of three unsigned integers.
inline unsigned int max3(unsigned int a,unsigned int b,unsigned int c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

/*
int cmpfunc (const void * a, const void * b){
	if ( ((edge*)a)->w>((edge*)b)->w ){
		return 1;
	}
	return -1;
}
*/

edgelist* readedgelist(char* input){
	unsigned e1=NLINKS;
	edgelist *g=malloc(sizeof(edgelist));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(input,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u\n", &(g->edges[g->e].s), &(g->edges[g->e].t))==2){//, &(g->edges[g->e].w))==3) {//Add one edge
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

	//qsort(g->edges,g->e,sizeof(edge),cmpfunc);
	return g;
}

// unionfind structure :
typedef struct {
	unsigned n;//number of objects
	unsigned *p;//parents
} unionfind;

unionfind* allocuf(unsigned n){
	unsigned i;
	unionfind* uf=malloc(sizeof(unionfind));
	uf->n=n;
	uf->p=malloc(n*sizeof(unsigned));
	for (i=0;i<n;i++){
		uf->p[i]=i;
	}
	return uf;
}

/*//create a singleton {x}
void MakeSet(unsigned x, unionfind *uf){
	uf->p[x]=x;
	uf->r[x]=0;
}
//Decided to do create directly all singletons
*/

//Merge the clusters of x and y returns 1 iff x and y belonged to the same cluster, 0 otherwise.
bool Union(unsigned x, unsigned y, unionfind *uf,omp_lock_t *lock){
	unsigned i, tmp;
	bool b;

	while (uf->p[x] != uf->p[y]){
		if (uf->p[x]<uf->p[y]){
			if (x==uf->p[x]){
    				omp_set_lock(&(lock[x]));
				b=0;
				if (x==uf->p[x]){
					uf->p[x]=uf->p[y];
					b=1;
				}
				omp_unset_lock(&(lock[x]));
				if (b)
					return 0;
			}
			tmp=uf->p[x];
			uf->p[x]=uf->p[y];
			x=tmp;
		}
		if (uf->p[x]>uf->p[y]){
			if (y==uf->p[y]){
    				omp_set_lock(&(lock[y]));
				b=0;
				if (y==uf->p[y]){
					uf->p[y]=uf->p[x];
					b=1;
				}
				omp_unset_lock(&(lock[y]));
				if (b)
					return 0;
			}
			tmp=uf->p[y];
			uf->p[y]=uf->p[x];
			y=tmp;
		}
	}
	return 1;
}

edgelist *alloctree(unsigned n){
	edgelist *el=malloc(sizeof(edgelist));
	el->edges=malloc((n-1)*sizeof(edge));
	el->n=0;
	el->e=0;
	return el;
}

edgelist* kruskal(edgelist* el){
	unsigned i,u,v;
	edgelist* elr=alloctree(el->n);
	unionfind *uf=allocuf(el->n);
	omp_lock_t *lock=malloc(uf->n*sizeof(omp_lock_t));

	for (i=0; i<uf->n; i++)
        	omp_init_lock(&(lock[i]));

	#pragma omp parallel for private(u,v) shared(uf,lock)
	for (i=0;i<el->e;i++){
		u=el->edges[i].s;
		v=el->edges[i].t;
		if (Union(u,v,uf,lock)==0){
    			#pragma omp atomic
			elr->e++;
			//elr->edges[elr->e++]=el->edges[i];
		}
	}
	return elr;
}

void printres(edgelist* el, char* output){
	FILE* file=fopen(output,"w");
	//double s=0;
	unsigned i;
	for (i=0;i<el->e;i++){
		fprintf(file,"%u %u\n",el->edges[i].s,el->edges[i].t);//,el->edges[i].w);
		//s+=el->edges[i].w;
	}
	fclose(file);
	//return s;
}


int main(int argc,char** argv){
	edgelist *el,*elr;

	omp_set_num_threads(atoi(argv[1]));

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edge list from file %s\n",argv[2]);
	el=readedgelist(argv[2]);


	printf("Number of nodes = %u\n",el->n);
	printf("Number of edges = %u\n",el->e);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing minimum spaning tree\n");

	elr=kruskal(el);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Printing result in file %s\n",argv[3]);

	//printres(elr,argv[3]);
	printf("Number of edges in resulting spaning tree: %u\n",elr->e);
	//printf("Sum of the weight in a minimum spaning tree: %le\n",s);

	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}
