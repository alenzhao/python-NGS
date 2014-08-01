
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
								int parent; char Aset;	float score; int coverage;
};




int mecscore(struct fragment* Flist,int fragments, char* h,float* ll, int* calls,int* miscalls)
{
								int j=0,k=0,good=0,bad=0,f=0; *ll =0; float p0,p1; *calls =0; *miscalls=0;
								for (f=0;f<fragments;f++)	
								{
																good=bad=0;p0=p1=0; //if (Flist[f].blocks ==1 && Flist[f].list[0].len ==1) continue;
																for (j=0;j<Flist[f].blocks;j++)
																{
																								*calls += Flist[f].list[j].len;
																								for (k=0;k<Flist[f].list[j].len;k++)
																								{
																																if (h[Flist[f].list[j].offset+k] == '-') continue;
																																if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good++; else																bad++;
																																if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) 
																																{ p0 += log10(1-Flist[f].list[j].pv[k]); p1 += log10(Flist[f].list[j].pv[k]); }
																																else
																																{ p1 += log10(1-Flist[f].list[j].pv[k]); p0 += log10(Flist[f].list[j].pv[k]); }
																								}
																}
																if (good < bad) *miscalls += good; else *miscalls += bad;
																if (p0 > p1) *ll += (p0 + log10(1+pow(10,p1-p0))); else *ll += (p1 + log10(1+pow(10,p0-p1))); 
																//*ll -= log10(2);
								}
}


int compute_fragscore(struct fragment* Flist,int f, char* h, float* ll)
{
								int j=0,k=0,good=0,bad=0;  float p0=0,p1=0;
								for (j=0;j<Flist[f].blocks;j++) 
								{
																for (k=0;k<Flist[f].list[j].len;k++) 
																{
																								if (h[Flist[f].list[j].offset+k] == '-') continue;
																								if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good++; else bad++;
																								if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) 
																								{ p0 += log10(1-Flist[f].list[j].pv[k]); p1 += log10(Flist[f].list[j].pv[k]); }
																								else
																								{ p1 += log10(1-Flist[f].list[j].pv[k]); p0 += log10(Flist[f].list[j].pv[k]); }
																}
								}
								if (p0 > p1) *ll= (p0 + log10(1+pow(10,p1-p0))); else *ll = (p1 + log10(1+pow(10,p0-p1)));
								//*ll -= log10(2);
								if (good < bad) return good; else return bad;
}

int update_fragscore(struct fragment* Flist,int f, char* h)
{
								int j=0,k=0,good=0,bad=0; float p0=0,p1=0;
								Flist[f].calls =0;
								for (j=0;j<Flist[f].blocks;j++) 
								{
																Flist[f].calls += Flist[f].list[j].len;
																for (k=0;k<Flist[f].list[j].len;k++) 
																{
																								if (h[Flist[f].list[j].offset+k] == '-') { fprintf(stdout,"fragment error"); exit(0);}
																								if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good++; else bad++;
																								if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) 
																								{ p0 += log10(1-Flist[f].list[j].pv[k]); p1 += log10(Flist[f].list[j].pv[k]); }
																								else
																								{ p1 += log10(1-Flist[f].list[j].pv[k]); p0 += log10(Flist[f].list[j].pv[k]); }
																}
								}
								if (p0 > p1) Flist[f].ll = (p0 + log10(1+pow(10,p1-p0))); else Flist[f].ll = (p1 + log10(1+pow(10,p0-p1)));
								//								fprintf(stdout,"frag %d score %f \n",f,Flist[f].ll);
								//Flist[f].ll -= log10(2);
								//								fprintf(stdout,"good %d bad %d \n",good,bad);
								if (good < bad) Flist[f].currscore = good; else Flist[f].currscore = bad;
}


void frag_cluster_initialize(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,char* h1,int snps,struct BLOCK* clist,int comps)
{
								int i=0,j=0,k=0,t=0,a,b,t1=0,t2=0,match=0, mismatch=0,diffs=0,flag=0,iter=0,alleles=0,assigned=0,first,last;
								char* thap = (char*)malloc(snps); 
								for (i=0;i<snps;i++) thap[i] = '-';
								for (i=0;i<snps;i++) h1[i] = '-';
								// cluster fragments within each connected component using pairwise overlap and initialize the starting haplotype from this solution 

								for (k=0;k<comps;k++)
								{
																for (i=0;i<clist[k].frags;i++) Flist[clist[k].flist[i]].clust = 'x'; // initialized 
																Flist[clist[k].flist[0]].clust = '0'; assigned=10; iter =0;
																while (assigned >0 || iter < 5)
																{
																								assigned =0; iter++;
																								for (i=0;i<clist[k].frags-1;i++)
																								{
																																t1 = clist[k].flist[i]; 	//if (Flist[t1].l == 'x') continue;
																																for (j=i+1;j<clist[k].frags;j++)
																																{
																																								last = Flist[t1].list[Flist[t1].blocks-1].offset+Flist[t1].list[Flist[t1].blocks-1].len; 
																																								t2 = clist[k].flist[j]; match =0; mismatch =0; 
																																								if (Flist[t2].list[0].offset >= last) break; 
																																								if (Flist[t1].clust != 'x' && Flist[t2].clust != 'x') continue;
																																								if (Flist[t1].clust == 'x' && Flist[t2].clust == 'x') continue;

																																								for (a=0;a<Flist[t1].blocks;a++)
																																								{
																																																for (b=0;b<Flist[t1].list[a].len;b++) thap[Flist[t1].list[a].offset+b] = Flist[t1].list[a].hap[b]; 
																																								}
																																								for (a=0;a<Flist[t2].blocks;a++)
																																								{
																																																for (b=0;b<Flist[t2].list[a].len;b++) { if (thap[Flist[t2].list[a].offset+b] == Flist[t2].list[a].hap[b]) thap[Flist[t2].list[a].offset+b] = 'M'; else thap[Flist[t2].list[a].offset+b] = 'm';}
																																								}
																																								for (a=0;a<Flist[t1].blocks;a++)
																																								{
																																																for (b=0;b<Flist[t1].list[a].len;b++) { if (thap[Flist[t1].list[a].offset+b] == 'M') match++; else if (thap[Flist[t1].list[a].offset+b] == 'm') mismatch++; }
																																								}
																																								if (match + mismatch ==0) continue;
																																								if (match != mismatch  && (Flist[t1].clust =='x' || Flist[t2].clust =='x')) assigned++;
																																								if (match - mismatch >= 1 && Flist[t1].clust != 'x') Flist[t2].clust = Flist[t1].clust;
																																								else if (mismatch - match >=1 && Flist[t1].clust =='0') Flist[t2].clust = '1';
																																								else if (mismatch -match >=1 && Flist[t1].clust =='1') Flist[t2].clust = '0';
																																								else if (match - mismatch >= 1 && Flist[t2].clust != 'x') Flist[t1].clust = Flist[t2].clust;
																																								else if (mismatch - match >=1 && Flist[t2].clust =='0') Flist[t1].clust = '1';
																																								else if (mismatch -match >=1 && Flist[t2].clust =='1') Flist[t1].clust = '0';
																																}
																								} //fprintf(stdout,"cluster %d frags %d assigned %d \n",k,clist[k].frags,assigned);

																}
								}
								for (i=0;i<fragments;i++) 
								{
																//														fprintf(stdout,"%d %d %c ",i,Flist[i].blocks,Flist[i].clust);for (t=0;t<Flist[i].blocks;t++) fprintf(stdout,"| %d %s ",Flist[i].list[t].offset,Flist[i].list[t].hap); fprintf(stdout,"\n");
																for (t1=0;t1<Flist[i].blocks;t1++)	
																{
																								for (t2=0;t2<Flist[i].list[t1].len;t2++)
																								{
																																if (Flist[i].clust == '0' && h1[Flist[i].list[t1].offset + t2] == '-') h1[Flist[i].list[t1].offset + t2] = Flist[i].list[t1].hap[t2]; 
																																if (Flist[i].clust == '1' && h1[Flist[i].list[t1].offset + t2] == '-' && h1[Flist[i].list[t1].offset + t2] =='1') h1[Flist[i].list[t1].offset + t2] ='0';
																																if (Flist[i].clust == '1' && h1[Flist[i].list[t1].offset + t2] == '-' && h1[Flist[i].list[t1].offset + t2] =='0') h1[Flist[i].list[t1].offset + t2] ='1';
																																// Flist[i].list[t1].pv[t2] = 0.05; bug here do not change pv HERE !!!!
																								}
																}
								}
								for (i=0;i<fragments;i++) { if (Flist[i].clust == 'x') Flist[i].clust = (char)(48+(int)(drand48()*2-0.00001)); }
								//for (i=0;i<snps;i++) { if (snpfrag[i].frags < 1 || h1[i] != '-') continue; h1[i] = '0';   sample_block(Flist,fragments,snpfrag,i,i,h1,snps); }
								for (i=0;i<snps;i++) { if (snpfrag[i].frags < 1 ) continue; h1[i] = '0';   sample_block(Flist,fragments,snpfrag,i,i,h1,snps); }


}


int sample_block(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int start, int end,char* h1,int snps)
{
								int j=0,t=0,t1=0,t2=0,flag=0,k=0,choice=0; char c;
								double parr0=0,parr1=0, sum=0,rd=0,p1,a,b;
								if (snpfrag[snpfrag[start].component].best_mec <= 0) return -1; if (snpfrag[start].frags < 1) return -1;
								for (j=snpfrag[start].ff;j<fragments;j++)
								{
																if (Flist[j].list[0].offset > end) break;
																if (Flist[j].component != snpfrag[start].component) continue;
																for (t1=0;t1<Flist[j].blocks;t1++)
																{
																								for (t2=0;t2<Flist[j].list[t1].len;t2++)
																								{
																																t = Flist[j].list[t1].offset + t2; 
																																if (t >= start && t <= end)
																																{
																																								b = log10(1-Flist[j].list[t1].pv[t2]) - log10(Flist[j].list[t1].pv[t2]);
																																								if (PMODEL ==0) 		b = log10(1-snpfrag[t].pv) - log10(snpfrag[t].pv);
																																								c = Flist[j].clust;  if (Flist[j].list[t1].hap[t2] != h1[t]) c=  (char)(49 - (int)(c-48));
																																								if (c =='0') parr0 += b; else parr0 -= b;
																																}
																								}
																}
								}
								if (parr0 >0 ) { parr1 -= parr0; parr0 = 0; } rd = drand48(); sum = pow(10,parr0) + pow(10,parr1);
								if (rd*sum >= pow(10,parr0))
								{
																for (j=start;j<=end;j++) { if (h1[j] == '-') continue; if (h1[j] =='0') h1[j] = '1'; else h1[j] = '0';  } 	choice = 1;
								}
								return choice;
}

void label_node(struct SNPfrags* snpfrag, int node,int comp)
{
								int i=0;
								if (snpfrag[node].component == -1) 
								{ 
																//											fprintf(stdout," called %d node edges %d %d \n",node,snpfrag[node].edges,comp);
																snpfrag[node].component = comp; snpfrag[comp].csize++; 
																for (i=0;i<snpfrag[node].edges;i++) label_node(snpfrag,snpfrag[node].elist[i].snp,comp); 
								}
}

void add_edges(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,int* components)
{
								int i=0,j=0,t=0,k=0,iter=0,maxdeg=0,avgdeg=0, mdelta=0,t1,t2,t3; char p,q;
								for (i=0;i<snps;i++) snpfrag[i].edges=0; 
								for (i=0;i<fragments;i++)
								{
																for (j=0;j<Flist[i].blocks;j++)
																{
																								for (k=0;k<Flist[i].list[j].len;k++) 
																								{
																																for (t=0;t<Flist[i].blocks;t++) 
																																{
																																								for (iter=0;iter<Flist[i].list[t].len;iter++)
																																								{
																																																if (Flist[i].list[j].offset+k == Flist[i].list[t].offset+iter) continue;
																																																if (Flist[i].list[j].offset +k - Flist[i].list[t].offset+iter > mdelta) mdelta = Flist[i].list[j].offset +k - Flist[i].list[t].offset+iter;
																																																snpfrag[Flist[i].list[t].offset+iter].elist[snpfrag[Flist[i].list[t].offset+iter].edges].snp = Flist[i].list[j].offset+k;
																																																snpfrag[Flist[i].list[j].offset+k].elist[snpfrag[Flist[i].list[j].offset+k].edges].frag = i;
																																																snpfrag[Flist[i].list[t].offset+iter].elist[snpfrag[Flist[i].list[t].offset+iter].edges].frag = i;
																																																snpfrag[Flist[i].list[t].offset+iter].elist[snpfrag[Flist[i].list[t].offset+iter].edges].p[0] = Flist[i].list[t].hap[iter];
																																																snpfrag[Flist[i].list[t].offset+iter].elist[snpfrag[Flist[i].list[t].offset+iter].edges].p[1] = Flist[i].list[j].hap[k] ;
																																																snpfrag[Flist[i].list[t].offset+iter].edges++; 
																																								}
																																}
																								}
																}
								}
								// sort all edges lists once for all by snp number 
								for (i=0;i<snps;i++)
								{
																for (j=0;j<snpfrag[i].edges-1;j++) 
																{
																								for (k=j+1; k<snpfrag[i].edges;k++)
																								{
																																if (snpfrag[i].elist[j].snp > snpfrag[i].elist[k].snp) 
																																{	
																																								t1 = snpfrag[i].elist[j].snp; snpfrag[i].elist[j].snp	= snpfrag[i].elist[k].snp; snpfrag[i].elist[k].snp = t1;
																																								t1 = snpfrag[i].elist[j].frag; snpfrag[i].elist[j].frag	= snpfrag[i].elist[k].frag; snpfrag[i].elist[k].frag = t1;
																																								p = snpfrag[i].elist[j].p[0]; snpfrag[i].elist[j].p[0]= snpfrag[i].elist[k].p[0]; snpfrag[i].elist[k].p[0] = p;
																																								p = snpfrag[i].elist[j].p[1]; snpfrag[i].elist[j].p[1]= snpfrag[i].elist[k].p[1]; snpfrag[i].elist[k].p[1] = p;
																																}
																								}
																}
								}
								for (i=0;i<snps;i++)
								{
																//								fprintf(stdout," snp %d edges %d || ",i,snpfrag[i].edges); for (j=0;j<snpfrag[i].edges;j++) fprintf(stdout,"%d ",snpfrag[i].elist[j]); fprintf(stdout,"\n"); getchar();
																if (snpfrag[i].edges > maxdeg) maxdeg = snpfrag[i].edges; avgdeg += snpfrag[i].edges;
																if (snpfrag[i].component != -1) continue; snpfrag[i].component = i;
																for (j=0;j<snpfrag[i].edges;j++) label_node(snpfrag,snpfrag[i].elist[j].snp,snpfrag[i].component);  
								}
								for (i=0;i<fragments;i++) Flist[i].component = snpfrag[Flist[i].list[0].offset].component; // each fragment has a component fixed 

								*components=0;	for (i=0;i<snps;i++) { if (snpfrag[i].component ==i && snpfrag[i].csize > 1) (*components)++; } fprintf(stdout," no of non-trivial connected components %d max-Degree %d avgdegree %f max delta %d\n",*components,maxdeg,(double)avgdeg/(double)snps,mdelta);        //getchar();
}

void update_snpfrags(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,int* components)
{
								int i=0,j=0,t=0,k=0,iter=0,maxdeg=0,avgdeg=0;
								for (i=0;i<snps;i++) { snpfrag[i].component = -1; snpfrag[i].csize=1; snpfrag[i].frags=0; snpfrag[i].edges=0; snpfrag[i].best_mec = 10000; snpfrag[i].best_mfr =0; snpfrag[i].no_best =0;} 
								for (i=0;i<fragments;i++)
								{
																for (j=0;j<Flist[i].blocks;j++)
																{
																								for (k=0;k<Flist[i].list[j].len;k++) 
																								{
																																snpfrag[Flist[i].list[j].offset+k].flist[snpfrag[Flist[i].list[j].offset+k].frags] = i;
																																snpfrag[Flist[i].list[j].offset+k].alist[snpfrag[Flist[i].list[j].offset+k].frags] = Flist[i].list[j].hap[k]; 
																																snpfrag[Flist[i].list[j].offset+k].frags++;
																																if (Flist[i].floppy ==1) continue;
																																for (t=0;t<Flist[i].blocks;t++) 
																																{
																																								for (iter=0;iter<Flist[i].list[t].len;iter++)
																																								{
																																																if (Flist[i].list[j].offset+k == Flist[i].list[t].offset+iter) continue;
																																																snpfrag[Flist[i].list[j].offset+k].edges++;
																																																snpfrag[Flist[i].list[t].offset+iter].edges++; 
																																								}
																																}
																								}
																}
								}
								fprintf(stdout,"updating snp data gibbs sampling procedure .......\n");
}

void output_current_solution(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,char* hap,char* best)
{
								// for each fragment, select if is going to be used in the solution or not 
								double score=0; int badfrags=0,components=0,i=0,j=0,k=0,t=0,s1,s2,b1,b2; char a,b; int pflag =0;
								for (i=0;i<fragments;i++) 
								{
																Flist[i].floppy =0; if (Flist[i].cons + Flist[i].bad ==0) { badfrags++; Flist[i].floppy = 1; continue; } 
																fprintf(stdout," %d len %d ratio %f \n",i,Flist[i].blocks,(double)Flist[i].cons/(double)(Flist[i].cons+Flist[i].bad));
																score = (double)Flist[i].cons/(double)(Flist[i].cons+Flist[i].bad);
																//		if (score < 0.9 && score > 0.1) { badfrags++; Flist[i].floppy = 1; }
								}
								fprintf(stderr,"printed fragment scores press enter to continue total %d bad %d \n\n",fragments,badfrags); //getchar();
								fprintf(stdout,"printed fragment scores press enter to continue total %d bad %d \n\n",fragments,badfrags); //getchar();
								//update_snpfrags(Flist,fragments,snpfrag,snps,&components); // connected components computed now 
								int happy=0, calls=0, miscalls=0, goodcalls=0,good,bad,comp,score_mfr=0,score_mec=0,frags_good=0; 
								for (i=0;i<snps;i++)	{ snpfrag[i].mfr =0; snpfrag[i].mec = 0; snpfrag[i].cfrags=0; snpfrag[i].mismatch=0;} 
								for (i=0;i<fragments;i++)
								{
																good =0; bad=0; t= calls; //if (Flist[i].floppy ==1) continue; 
																comp = Flist[i].component; if (comp <0) continue; snpfrag[comp].cfrags++;
																for (j=0;j<Flist[i].blocks;j++) 
																{
																								calls += Flist[i].list[j].len;    snpfrag[comp].mec += Flist[i].list[j].len;
																								for (k=0;k<Flist[i].list[j].len;k++) 
																								{
																																if ( hap[Flist[i].list[j].offset+k] == '-') continue;
																																if (Flist[i].list[j].hap[k] == best[Flist[i].list[j].offset+k]) good++; else bad++;
																								}
																} 
																if (good < bad && good > 0 && bad > 0 ) snpfrag[comp].mfr += good; else if (good > 0 && bad > 0) snpfrag[comp].mfr += bad;
								} 
								for (i=0;i<snps;i++)	
								{
																if (snpfrag[i].csize <= 1) continue; j=snpfrag[i].csize; t=i;
																while (j > 0 && t < snps) {	if (snpfrag[t].component == snpfrag[i].component) j--; t++; }
																fprintf(stderr,"BLOCK: offset %d len %d phased %d MFR %d frags %d MEC %d calls %d pairs of mismatching fragments %d\n",i+1,t-i,snpfrag[i].csize,snpfrag[i].best_mfr,snpfrag[i].cfrags,snpfrag[i].best_mec,snpfrag[i].mec,snpfrag[i].mfr); 
																fprintf(stdout,"BLOCK: offset %d len %d phased %d MFR %d frags %d MEC %d calls %d pairs of mismatching fragments %d\n",i+1,t-i,snpfrag[i].csize,snpfrag[i].best_mfr,snpfrag[i].cfrags,snpfrag[i].best_mec,snpfrag[i].mec,snpfrag[i].mfr); 
																if (snpfrag[i].best_mec > snpfrag[i].mfr) { for (j=i;j<t;j++) fprintf(stdout,"%c",hap[j]); fprintf(stdout,"\n"); } 
																if (snpfrag[i].best_mec > snpfrag[i].mfr) { for (j=i;j<t;j++) fprintf(stdout,"%c",best[j]); fprintf(stdout,"\n"); } 
																if (pflag) fprintf(stdout,"AARON BLOCK: offset: %d len: %d phased: %d \n",i+1,t-i,snpfrag[i].csize); 
																for (j=i;j<t;j++) 
																{
																								if (hap[j]=='0') { a ='0'; b= '1';} else if (hap[j] =='1') { a = '1'; b = '0';} else { a='-'; b =  '-'; } 
																								if (pflag) { if (snpfrag[i].component != snpfrag[j].component) fprintf(stdout,"AARON %5d -  - \n",j+1); else fprintf(stdout,"AARON %5d %c  %c \n",j+1,a,b);}
																} if (pflag) fprintf(stdout,"AARON ********\n");
																score_mfr += snpfrag[i].best_mfr; score_mec += snpfrag[i].best_mec;
																fprintf(stdout,"\n**************************************************\n");
								}

								fprintf(stderr," outputted current solution with %d components....... frags %d MFR %d MEC %d calls \n\n",components,frags_good,score_mfr,score_mec); //getchar();	
								fprintf(stdout," outputted current solution with %d components....... frags %d MFR %d MEC %d calls \n\n",components,frags_good,score_mfr,score_mec); //getchar();	


}
int print_hapfile(struct BLOCK* clist,int blocks,char* h1,struct fragment* Flist,int fragments,struct SNPfrags* snpfrag, char* fname,int score)
{
								// print a new file containing one block phasing and the corresponding fragments 
								int i=0,j=0,t=0,k=0; char c;
								char fn[200]; sprintf(fn,"%s-%d.phase",fname,score); 
								FILE* fp; 		 fp = fopen(fn,"w");

								for (i=0;i<blocks;i++)
								{
																fprintf(fp,"BLOCK: offset: %d len: %d phased: %d \n",clist[i].offset+1,clist[i].length,clist[i].phased); 	
																for(k=0;k<clist[i].length;k++) 
																{
																								t=clist[i].offset+k; 
																								if (clist[i].haplotype[k] =='-') fprintf(fp,"%s\t-\t- \n",snpfrag[t].id);
																								else
																								{
																																if (h1[t] =='0') c= '1'; else if (h1[t] =='1') c = '0'; 
																																fprintf(fp,"%s\t%c\t%c \n",snpfrag[t].id,h1[t],c); 
																								}
																}
																if (i < blocks-1) fprintf(fp,"******** \n");
								}
								fclose(fp);
}


int print_block_frags(struct BLOCK* blist, int block,char* aaron,char* h1,int* bn,struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,char* fragfile)
{
								// print a new file containing one block phasing and the corresponding fragments 
								int i=0,j=block,pflag=1,offset = blist[j].offset,k=0,t=0; char c;
								fprintf(stdout," starting to print block %d length %d phased %d offset %d\n",block,blist[block].length,blist[block].phased,offset);

								FILE* fp;
								//char fname[100]; sprintf(fname,"%s-BLOCK-%d-offset_%d.phase",fragfile,block,offset); 
								//fp = fopen(fname,"w");
								fp = stdout;
								if (pflag) fprintf(fp,"BLOCK: offset: %d len: %d phased: %d \n",1,blist[j].length,blist[j].phased); 	
								for(k=0;k<blist[j].length;k++) 
								{
																t=blist[j].offset+k; 
																if (bn[t] != bn[blist[j].offset] && pflag ==1) fprintf(fp,"%s\t-\t- \n",snpfrag[t].id);
																else
																{
																								if (h1[t] =='-') c = '-'; else if (h1[t] =='0') c= '1'; else if (h1[t] =='1') c = '0'; 
																								if (pflag) fprintf(fp,"%s\t%c\t%c \n",snpfrag[t].id,h1[t],c); 
																}
								}
								//								if (pflag) fprintf(fp,"******** \n");
								//fclose(fp); sprintf(fname,"%s-BLOCK-%d-offset_%d.matrix",fragfile,block,offset);  fp = fopen(fname,"w");

								int frags=0,start=-1;
								for (j=0;j<fragments;j++)
								{
																if (Flist[j].list[0].offset < blist[block].offset ) continue; 
																if (start == -1) start = j;
																if (Flist[j].component != snpfrag[offset].component) continue;
																if (Flist[j].list[0].offset >= offset + blist[block].length) break;
																frags++;			
								}								
								fprintf(fp,"%d %d \n",frags+1,blist[block].length);
								for (j=start;j<fragments;j++)
								{
																if (Flist[j].component != snpfrag[offset].component) continue;
																if (Flist[j].list[0].offset >= offset + blist[block].length) break;
																fprintf(fp,"%s ",Flist[j].id); for (k=0;k<Flist[j].blocks;k++) fprintf(fp,"%d %s ",Flist[j].list[k].offset-offset+1,Flist[j].list[k].hap); fprintf(fp,"\n");
								} 
								//fclose(fp);

								//							fprintf(stderr,"printing bad block \n"); exit(0);

}

int print_blocks(struct BLOCK* blist, int blocks,char* aaron,char* h1,char* current,int* bn,struct fragment* Flist, int fragments, int mecscore,struct SNPfrags* snpfrag) // h1 is the current best, aaron is initial and current is current sampled haplotype 
{
								fprintf(stdout,"START OF OUTPUT OF HAPLOTYPE PHASED BLOCKS \n\n");
								int j=0; for (j=0;j<blocks;j++) print_block(blist,j,aaron,h1,current,bn,Flist,fragments,mecscore,snpfrag);
}

int print_block(struct BLOCK* blist, int block,char* aaron,char* h1,char* current,int* bn,struct fragment* Flist, int fragments, int mecscore,struct SNPfrags* snpfrag) // h1 is the current best, aaron is initial and current is current sampled haplotype 
{
								int i=0,j=0,k=0,l=0, w=0,am =0, nm=0,cmm=0,amm=0,nmm=0,common=0,t=0,good=0,bad=0,AM=0,GM=0,cm=0,a=0,diffs=0,op=0,delta=0,pos=0,mec1,mec2,mec1sum=0,mec2sum=0,flag=0,pflag=0,window = 115,first,last,frags=0,omatch=0,nmatch=0,lf=0,islands=0; 
								char ht[200]; char c;float ll=0;
								for(i=0;i<fragments;i++) Flist[i].floppy =0;
								j=block;
								{
																diffs =0; for (i=blist[j].offset;i<blist[j].offset+blist[j].length;i++) { if (h1[i] != aaron[i] && bn[i]==bn[blist[j].offset]) diffs++; }
																pos = 0; for (i=blist[j].offset;i<blist[j].offset+blist[j].length;i++) { if (bn[i]==bn[blist[j].offset]) pos++; }
																islands = 0; 
																//for (i=blist[j].offset;i<blist[j].offset+blist[j].length;i++) { if (bn[i]==bn[blist[j].offset] && snpfrag[i].FIRST == i) islands++; }
																frags=0;omatch=0;nmatch=0;mec1sum=0;mec2sum=0;
																for(l=0;l<blist[j].frags;l++) 
																{
																								i = blist[j].flist[l]; 
																								mec1 = compute_fragscore(Flist,i,aaron,&ll); 	mec1sum += mec1; Flist[i].sc1 = mec1; if (mec1==0) omatch++;
																								mec2 = compute_fragscore(Flist,i,h1,&ll); 	mec2sum += mec2; Flist[i].sc2 = mec2; if (mec2==0) nmatch++;
																}

																for (a=0;a<130;a++) fprintf(stdout,"=");
																fprintf(stdout,"\nBLOCK: OFFSET %9d length %5d phased %5d \n",blist[j].offset+1,blist[j].length,blist[j].phased);
																if (mec2sum -mec1sum >= 1) { print_block_frags(blist,j,aaron,h1,bn,Flist,fragments,snpfrag,"S");	fprintf(stdout,"\nBLOCK: OFFSET %9d length %5d phased %5d \n",blist[j].offset+1,blist[j].length,blist[j].phased); } 
																fprintf(stdout,"frags %d ORIG-match %d new-match %d ORIG-MEC %d NEW-MEC %d ISLANDS in block %d \n",blist[j].frags,omatch,nmatch,mec1sum,mec2sum,islands);
																if (mec2sum < mec1sum) fprintf(stdout,"NEW SOLUTION IS BETTER THAN ORIGINAL.... \n"); 
																if (mec2sum > mec1sum) fprintf(stdout,"NEW SOLUTION IS WORSE %d THAN ORIGINAL.... \n",mec2sum-mec1sum); 
																if (mec2sum ==0) return 1;
																if ( ( diffs == 0 || pos == diffs) ) return 1;

																first = blist[j].offset; last = first+ window;
																while (1)
																{
																								if (last >= blist[j].offset + blist[j].length) last = blist[j].offset + blist[j].length;
																								for (a=0;a<130;a++) fprintf(stdout,"-"); fprintf(stdout,"\n\n");
																								for (i=first;i<last;i++)
																								{
																																if (bn[i]==bn[blist[j].offset] && snpfrag[i].frags > 9) fprintf(stdout,"+");
																																else if (bn[i]==bn[blist[j].offset]) fprintf(stdout,"%d",snpfrag[i].frags);
																																else fprintf(stdout,"-");
																								}        fprintf(stdout,"  frag coverage \n");
																								for (i=first;i<last;i++)
																								{
																																if (bn[i]==bn[blist[j].offset]) fprintf(stdout,"%c",snpfrag[i].island); else fprintf(stdout,"-");
																								}        fprintf(stdout,"  island info \n");
																								for (i=first;i<last;i++)
																								{
																																if (h1[i] != aaron[i] && bn[i]==bn[blist[j].offset]) fprintf(stdout,"x"); else fprintf(stdout,".");
																								}fprintf(stdout,"\n");


																								for (i=first;i<last;i++)
																								{
																																if (bn[i]==bn[blist[j].offset]) fprintf(stdout,"%c",current[i]); else fprintf(stdout,"-");
																								}
																								fprintf(stdout,"  SAMPLE HAP.\n");
																								for (i=first;i<last;i++)
																								{
																																if (current[i] != aaron[i] && bn[i]==bn[blist[j].offset]) fprintf(stdout,"x"); else fprintf(stdout,".");
																								}fprintf(stdout,"\n");
																								for (i=first;i<last;i++)
																								{
																																if (bn[i]==bn[blist[j].offset]) fprintf(stdout,"%c",aaron[i]); else fprintf(stdout,"-");
																								}
																								fprintf(stdout,"  INITIAL HAP.\n");
																								for (i=first;i<last;i++)
																								{
																																if (h1[i] != aaron[i] && bn[i]==bn[blist[j].offset]) fprintf(stdout,"x"); else fprintf(stdout,".");
																								}fprintf(stdout,"\n");
																								for (i=first;i<last;i++)
																								{
																																if (bn[i]==bn[blist[j].offset]) fprintf(stdout,"%c",h1[i]); else fprintf(stdout,"-");
																								} fprintf(stdout,"  CUR-BEST HAP.\n\n");

																								for (i=snpfrag[first].ff;i<fragments;i++)
																								{
																																if (Flist[i].list[0].offset >= last) break;
																																if (Flist[i].component != snpfrag[blist[j].offset].component) continue;
																																if (Flist[i].sc1 ==0 &&  Flist[i].sc2 ==0) continue;
																																for (k=0;k<window;k++) ht[k] = '-';
																																for (k=0;k<Flist[i].blocks;k++)
																																{
																																								for (t=0;t<Flist[i].list[k].len;t++)
																																								{
																																																if (Flist[i].list[k].offset + t >= first && Flist[i].list[k].offset + t < last) ht[Flist[i].list[k].offset + t-first] = Flist[i].list[k].hap[t];
																																								}
																																} Flist[i].floppy = 1;
																																for (k=0;k<window && k < last-first;k++) fprintf(stdout,"%c",ht[k]);   fprintf(stdout,"f %d %d %d \n",i,Flist[i].sc1,Flist[i].sc2);
																								}fprintf(stdout,"\n");
																								while (1)  
																								{
																																good =0;
																																for (k=0;k<window;k++) ht[k] = '.'; lf =-10;
																																for (i=snpfrag[first].ff;i<fragments;i++)
																																{
																																								if (Flist[i].list[0].offset >= last) break;
																																								if (Flist[i].component != snpfrag[blist[j].offset].component) continue;
																																								if (Flist[i].floppy ==1) continue;
																																								if (Flist[i].list[0].offset-first <= lf+2 && lf > 0) continue; 
																																								if (lf >0) ht[lf+1] = ')';
																																								for (k=0;k<Flist[i].blocks;k++)
																																								{
																																																for (t=0;t<Flist[i].list[k].len;t++)
																																																{
																																																								if (Flist[i].list[k].offset + t >= first && Flist[i].list[k].offset + t < last ) ht[Flist[i].list[k].offset + t-first] = Flist[i].list[k].hap[t];
																																																								if (Flist[i].list[k].offset + t >= first && Flist[i].list[k].offset + t < last) lf = Flist[i].list[k].offset + t-first; 
																																																}
																																								} Flist[i].floppy = 1; good++;
																																} 
																																if (lf > 0) ht[lf+1] = ')';
																																flag =0;
																																for (k=0;k<window && k < last-first;k++) 
																																{
																																								if (ht[k] == '1' || ht[k] =='0') flag = 1;
																																								else if (ht[k] == '.' && flag ==1) ht[k] = '=';
																																								else if (ht[k] == ')' ) flag =0;
																																}
																																t=0; for (k=0;k<window && k < last-first;k++) { if (ht[k] == ')') ht[k] = '.';  if (ht[k] != '.') t++;}
																																if (t>0) { for (k=0;k<window && k < last-first;k++)  fprintf(stdout,"%c",ht[k]);   fprintf(stdout," perfect match\n"); } 
																																if (good ==0) break;
																								}

																								first += window; last = first + window;
																								if (first >= blist[j].offset + blist[j].length) break;

																}

								}
}

