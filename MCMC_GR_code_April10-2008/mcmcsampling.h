
/************************************* CODE FOR DOINF TREE BASED EXACT SAMPLING **********************************************/

int switcherror_estimates(struct SNPfrags* snpfrag, int* H0011, int* H0110,char* h1, struct BLOCK* clist, int components,int Z);

int update_posterior(struct fragment* Flist,int fragments,char* h,int Z);

int sample_pair_hap(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1,int compute);

int sample_local_hap(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1,char** HAP,float* prob,float* probtemp,int compute);

int compute_tree_prob(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp);

int update_tree_prob(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp);

int sample_tree(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp,int ML); 

int construct_tree(struct BLOCK* clist,int k,int parent, int st);  // PREORDER TRAVERSAL ---> BINARY TREE 


int sample_haplotype(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k,char* h1,char** HAP,float* prob,float* probtemp);

int sample_block(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int start, int end,char* h1,int snps);


int compute_flist(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node, int mode);
int NJtree(struct SNPfrags* snpfrag, int snps,struct BLOCK* clist, int block,FILE* mcout,int* nodes);
 
int compute_updatelist(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node);
								// find all other subsets (nodes) that are affected by flipping of this subset 
int compute_flipprob(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1,int update);

int sample_hap_rest(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1);

int partitions(struct SNPfrags* snpfrag, int snps,int* slist, int N,int min,struct fragment* Flist,int fragments,FILE* mcout,int* nodes,char* hap);
int min_cut(struct SNPfrags* snpfrag, int snps,int* slist, int N,char* hap);


								// compute min-cut of graph represented by vertices in slist of size N 
								// this function returns a list of SUBSETS in order that represent consecutive MIN-CUTS
								// struct edge* telist; int tedges; // temporary edge list and number of edges for MIN CUT computation 
