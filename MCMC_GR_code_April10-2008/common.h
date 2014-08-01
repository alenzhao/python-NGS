
//Tue May 29 23:13:29 PDT 2007
//


extern int GIBBS1;
extern int TREE;
extern int AARON;
extern int STDERR;
extern int RANDOM_START;
extern int USE_HAP;
extern int ERROR;
extern int QV;
extern int MINCUT;
extern int burnin;
extern int thinrate;
extern int MCMCruns;
extern int SINGLE;
extern int UPDATE;
extern int UPDATETREE;
extern int EXAMPLE; extern int ne; extern int MF;
extern int PMODEL;







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

