
//Tue May 29 23:13:29 PDT 2007
//

//#include "common.h"

int compare_haps(struct BLOCK* clist, int components, char* orig, char* h1,struct SNPfrags* snpfrag, int snps);
int compare_Flist_hap(struct SNPfrags* snpfrag, int snps, struct fragment* Flist,int fragments,char* h,int Z,int QV);
int mutate_Flist(struct fragment* Flist,int fragments,double errrate);
int correct_fragment(struct fragment* Flist,int f, char* h);
int mecscore(struct fragment* Flist,int fragments, char* h,float* ll, int* calls,int* miscalls);
int compute_fragscore(struct fragment* Flist,int f, char* h, float* ll);
int update_fragscore(struct fragment* Flist,int f, char* h);
void frag_cluster_initialize(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,char* h1,int snps,struct BLOCK* clist,int comps);

void label_node(struct SNPfrags* snpfrag, int node,int comp);
void add_edges(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,int* components);
void update_snpfrags(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,int* components);
void output_current_solution(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,char* hap,char* best);

int print_hapfile(struct BLOCK* clist,int blocks,char* h1,struct fragment* Flist,int fragments,struct SNPfrags* snpfrag, char* fname,int score);
int print_block_frags(struct BLOCK* blist, int block,char* aaron,char* h1,int* bn,struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,char* fragfile);

int print_blocks(struct BLOCK* blist, int blocks,char* aaron,char* h1,char* current,int* bn,struct fragment* Flist, int fragments, int mecscore,struct SNPfrags* snpfrag);

int print_block(struct BLOCK* blist, int block,char* aaron,char* h1,char* current,int* bn,struct fragment* Flist, int fragments, int mecscore,struct SNPfrags* snpfrag);

void generate_example_2(int n);

void generate_example(int n,int M);

int phase_score(char* hap,int i, int j, char* p );






