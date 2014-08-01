
//Tue May 29 23:13:29 PDT 2007
//
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

int PMODEL = 1;
struct NODE // node of tree 
{
								int offset; short length; short phased; char* hap; char* bitvec; int frags; int* flist; // for siblings 
								int ch1,ch2,parent; char leaf;
								float prob,p01,p00;short H;
								float P00,P01; // temporary 
								//char** HAP; float* pdf; // prob is probability of current node haplotype, HAP & pdf are for leaves
								int* sflist; int sfrags; int flips,potflips;
								int* updatelist; // list of other SNP subsets that are affected by flipping this subset and hence should be updated 
								int updatenodes;
};

#define MAXBSIZE 100
struct block
{
								int offset; char* hap; short len; float* pv; char* qv; float* post; 
};

struct BLOCK
{
								int offset; int length; int phased; char* haplotype; int* flist; int frags; int islands;
								struct NODE* tree; int nodes; // nodes is the number of nodes in tree  
								int MEC, bestMEC,lastMEC,treeMEC,treebestMEC; 
								float LL, bestLL,treeLL,treebestLL;  // log likellihood scores 
								int calls; char treecompute,dealloc;
};

struct fragment
{
								char id[20]; short blocks; struct block* list; char clust; // l indicates if it belongs to the haplotype or it's complement
								int cons,bad; int floppy; // floppy indicates a bad fragment
								char single;
								int component; int sc1,sc2; int currscore,calls; float ll;
};

struct edge
{
								int snp; int frag; char p[2]; float w; 
};

struct SNPfrags 
{
								char id[20]; int* flist; int frags; char* alist; // alist is the alleles corresponding to each fragment that covers this SNP
								int component; int edges; // those that span the interval between this snp and the next snp 
								int csize; int cfrags, mfr, mec; int blockno; int best_mec; int best_mfr; int mismatch; int no_best; int calls; short all_best;
								struct edge* elist;  int FIRST; int LAST; int ff; double pv; char island; double node_id;
								int bcomp; // index of clist to which this snp belongs: reverse mapping  
								struct edge* telist; int tedges; // temporary edge list and number of edges for MIN CUT computation 
								// FIRST is the first snp that belongs to an island and LAST is the last one, in between there could be gaps 
								int parent; char Aset;	float score;
};

/************************** FUNCTIONS FOR MANIPULATING SNPFRAG etc *****************************/

void label_node(struct SNPfrags* snpfrag, int node,int comp);

void add_edges(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,int* components);

void update_snpfrags(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,int* components);

void output_current_solution(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,char* hap,char* best);

int print_hapfile(struct BLOCK* clist,int blocks,char* h1,struct fragment* Flist,int fragments,struct SNPfrags* snpfrag, char* fname,int score);
								// print a new file containing one block phasing and the corresponding fragments 

int print_block_frags(struct BLOCK* blist, int block,char* aaron,char* h1,int* bn,struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,char* fragfile);

int print_blocks(struct BLOCK* blist, int blocks,char* aaron,char* h1,char* current,int* bn,struct fragment* Flist, int fragments, int mecscore,struct SNPfrags* snpfrag);
 // h1 is the current best, aaron is initial and current is current sampled haplotype 

int print_block(struct BLOCK* blist, int block,char* aaron,char* h1,char* current,int* bn,struct fragment* Flist, int fragments, int mecscore,struct SNPfrags* snpfrag);
 // h1 is the current best, aaron is initial and current is current sampled haplotype 

int compare_haps(struct BLOCK* clist, int components, char* orig, char* h1,struct SNPfrags* snpfrag, int snps)
int compare_Flist_hap(struct SNPfrags* snpfrag, int snps, struct fragment* Flist,int fragments,char* h,int Z,int QV)

int mutate_Flist(struct fragment* Flist,int fragments,double errrate)
int correct_fragment(struct fragment* Flist,int f, char* h)
int mecscore(struct fragment* Flist,int fragments, char* h,float* ll, int* calls,int* miscalls)
int compute_fragscore(struct fragment* Flist,int f, char* h, float* ll)
int update_fragscore(struct fragment* Flist,int f, char* h)



/*************************** CODE FOR DOING MCMC  SAMPLING **********************************************/

int update_posterior(struct fragment* Flist,int fragments,char* h,int Z)
int gibbs_sampling(char* frags,char* sol,int fraction);
int sample_block(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int start, int end,char* h1,int snps);


int sample_pair_hap(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1,int compute)
int sample_local_hap(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1,char** HAP,float* prob,float* probtemp,int compute)

int compute_tree_prob(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp)

int update_tree_prob(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp)

int sample_tree(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp,int ML) // ML ==2 iimplies maximum likelihood else sample 
								// sample haplotype for block k using the tree clist[k].tree 
								//fprintf(stdout,"node %d offset %d phased %d child1 %d child2 %d frags %d\n",node,clist[k].tree[node].offset,clist[k].tree[node].phased,clist[k].tree[node].ch1,clist[k].tree[node].ch2,clist[k].tree[node].frags);

int construct_tree(struct BLOCK* clist,int k,int parent, int st);
 // PREORDER TRAVERSAL ---> BINARY TREE 


int switcherror_estimates(struct SNPfrags* snpfrag, int* H0011, int* H0110,char* h1, struct BLOCK* clist, int components,int Z);

int sample_haplotype(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k,char* h1,char** HAP,float* prob,float* probtemp);


void frag_cluster_initialize(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,char* h1,int snps,struct BLOCK* clist,int comps);

int sample_block(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int start, int end,char* h1,int snps);


