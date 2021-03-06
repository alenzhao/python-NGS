// modified April 9 2008 
//
//
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#include "common.h"

#include "fragmatrix.h"
#include "mcmcsampling.h"


int GIBBS =1;
int TREE =0;
int AARON =1;
int STDERR =1;
int RANDOM_START=0;
int USE_HAP = 1;
int ERROR = 1;
int QV = -1;
int MINCUT =1;  // use tree decomposition or Gamma_1 (0) 
int burnin = 20000;
int thinrate = 100;
int MCMCruns=0;
int SINGLE =0; // gamma_1 
int UPDATE =1;
int UPDATETREE =1;
int EXAMPLE = 0; int ne = 20; int MF = 3;
int PMODEL =1;



int tree_sampling(char* frags,char* sol,int fraction);

int main(int argc, char** argv)
{
		      printf("HASH program: MCMC sampling based Haplotype assembly.... \n");
								char frags[100]; char sol[100]; 
								if (EXAMPLE) 
								{
															 UPDATE =0; UPDATETREE=0;  
																ne = atoi(argv[1]); if (argc > 2) MF = atoi(argv[2]); 
																generate_example(ne,MF);
																sprintf(frags,"example-%d.matrix",ne); sprintf(sol,"example-%d.phase",ne);
																fprintf(stdout,"example file %s \n",frags); 
																tree_sampling(frags,sol,0);
																return 1;
								} 

								int flip =0;			if (argc >= 4) flip = atoi(argv[3]);
								if (argc <3) return 1;
								tree_sampling(argv[1],argv[2],flip); 
}

int tree_sampling(char* frags,char* sol,int fraction)
{
								char newfrags[200]; char qvfile[200];
								char command[500]; char command1[500];
								sprintf(command,"more %s | awk 'BEGIN {rows=1; snps=0;} { if ($3 == \"\") snps = $2; else if ($4 != NULL || length($3) > 1) {rows =rows+1; print $0;} } END { print rows,snps;} ' |  awk ' { flag =0; f =4; for (i =3; i<= 100; i+=1) {  if ($i == NULL && flag ==0  && $i != 0) { f= i; flag =1; }}  t = f-2; if (f >=4) print t/2,$0;  else print $0; }' | sort -g -k 3 -k 1 > %s.SORTED ;",frags,frags);
								//        sprintf(command," awk ' { flag =0; f =4; for (i =3; i<= 100; i+=1) {  if ($i == NULL && flag ==0  && $i != 0) { f= i; flag =1; }}  t = f-2; if (f >=4) print t/2,$0;  else print $0; }' %s | sort -g -k 3 -k 1 > %s.SORTED ;",frags,frags);
								sprintf(qvfile,"%s_qv",frags); fprintf(stderr,"qv file %s \n",qvfile); 
								FILE* ft = fopen(qvfile,"r"); if (ft == NULL) QV = -1; else fclose(ft); 

								if (frags[strlen(frags)-1] == 'd') { fprintf(stdout," filename has to end in SORTED or .matrix...... exittting \n");exit(0);}
								if (frags[strlen(frags)-1] == 'x')
								{
																fprintf(stdout,"executing awk command to sort fragment file  %s \n",frags); 
																system(command); sprintf(newfrags,"%s.SORTED",frags);
																if (QV != -1) 
																{
																								sprintf(command,"more %s_qv | awk 'BEGIN {rows=1; snps=0;} { if ($3 == \"\") snps = $2; else if ($4 != NULL || length($3) > 1) {rows =rows+1; print $0;} } END { print rows,snps;} ' |  awk ' { flag =0; f =4; for (i =3; i<= 100; i+=1) {  if ($i == NULL && flag ==0  && $i != 0) { f= i; flag =1; }}  t = f-2; if (f >=4) print t/2,$0;  else print $0; }' | sort -g -k 3 -k 1 > %s.SORTED_qv ;",frags,frags);
																								system(command); 
																								sprintf(qvfile,"%s.SORTED_qv",frags); fprintf(stderr,"qv file %s \n",qvfile); 
																								fprintf(stderr,"qv file %s \n",qvfile); //getchar();
																}
								}
								else sprintf(newfrags,"%s",frags);
								// IMP NOTE: all SNPs start from 1 instead of 0 and all offsets are 1+
								fprintf(stdout,"calling gibbs sampling procedure\n");
								int fragments=0,snps=0,len,phased=0,iter=0,components=0,pflag=0,t1=0,t2=0,flag=0;
								int i=0,j=0,k=0,t,l,biter=0,offset,blocks,type=0,component;
								char c1,c2; 
								int* slist; int min=1; FILE* MCout; int maxphased = 20000;  char tchar[20];
								int	haps = pow(2,min); 	char** HAP = (char**)malloc(sizeof(char*)*haps);
								for (i=0;i<haps;i++){ HAP[i] = (char*)malloc(min+1);	j=i; for (t=0;t<min;t++) {if (j%2 ==0) HAP[i][t] = '0'; else HAP[i][t] ='1'; j = j/2; } HAP[i][min] = '\0';} // 000 100 010 110  001  101  011  111  reverse order 
								float* prob = (float*)malloc(sizeof(float)*haps); float q;
								float* ptemp0 = (float*)malloc(sizeof(float)*haps); 
								float* ptemp1 = (float*)malloc(sizeof(float)*haps); 
								int mechist[50]; for (i=0;i<50;i++) mechist[i] =0; int flip,prevscore,allflips=0,realblockflips=0,blockflips=0,realallflips=0;


								FILE* ff = fopen(newfrags,"r"); if (ff == NULL) { fprintf(stderr,"couldn't open fragment file \n"); exit(0);}
								fscanf(ff,"%d %d \n",&fragments,&snps);  fragments--; // fragments is actually the number lines 
								fprintf(stdout,"fragments %d snps %d \n",fragments,snps); 
								fclose(ff);

								/****************************** READ FRAGMENT MATRIX*************************************************/
								ff = fopen(newfrags,"r"); if (ff == NULL) { fprintf(stderr,"couldn't open fragment file \n"); exit(0);}
								char buffer[1024]; char id[20]; char block[200]; char ch;
								struct fragment* Flist = (struct fragment*)malloc(sizeof(struct fragment)*fragments); 
								ch = fgetc(ff); while (ch != '\n') ch= fgetc(ff); 
								for (i=0;i<fragments;i++) 
								{
																j=0; ch = fgetc(ff); while (ch !='\n') { buffer[j] = ch; j++; ch = fgetc(ff); } buffer[j] = '\0';
																k=0; t=0; type=0; while (k < j)
																{
																								while (buffer[k] !=' ' && k < j && buffer[k] !='\0') { block[t] = buffer[k];t++; k++; } k++; 
																								while (buffer[k] ==' ' && k < j) k++; block[t] = '\0'; 
																								if (type ==0) 
																								{
																																blocks =0; for (l=0;l<t;l++) { blocks = 10*blocks + (int)(block[l]-48); } type =1; Flist[i].blocks = blocks;
																																Flist[i].list = (struct block*)malloc(sizeof(struct block)*(blocks));    biter=0;
																								}
																								else if (type ==1)
																								{
																																strcpy(Flist[i].id,block); type =2;
																								}
																								else if (type ==2)
																								{
																																offset = 0; for (l=0;l<t;l++) { offset = 10*offset + (int)(block[l]-48); } type =3; Flist[i].list[biter].offset = offset-1;
																								}
																								else if (type ==3)
																								{
																																Flist[i].list[biter].hap = (char*)malloc(t+1); strcpy(Flist[i].list[biter].hap,block); Flist[i].list[biter].len = t;
																																Flist[i].list[biter].pv = (float*)malloc(sizeof(float)*t);
																																for (t1=0;t1<Flist[i].list[biter].len;t1++) Flist[i].list[biter].pv[t1] = 0.01;
																																Flist[i].list[biter].qv = (char*)malloc(t+1);
																																for (t1=0;t1<Flist[i].list[biter].len;t1++) Flist[i].list[biter].qv[t1] = '0';
																																Flist[i].list[biter].post = (float*)malloc(sizeof(float)*t); // how many times it matches 
																																for (t1=0;t1<Flist[i].list[biter].len;t1++) Flist[i].list[biter].post[t1] =0;
																																type =2; biter++;
																								}
																								t=0;
																} 
								} fclose(ff);
								/****************************** READ FRAGMENT MATRIX*************************************************/

								/****************************** READ FRAGMENT QUALITY FILE*************************************************/
								ff = fopen(qvfile,"r"); 
								if (ff == NULL || QV == -1) fprintf(stderr,"couldn't open fragment QV file \n");
								else
								{
																fprintf(stderr,"reading fragment quality file \n");
																ch = fgetc(ff); while (ch != '\n') ch= fgetc(ff); 
																for (i=0;i<fragments;i++) 
																{
																								biter=0; j=0; ch = fgetc(ff); while (ch !='\n') { buffer[j] = ch; j++; ch = fgetc(ff); } buffer[j] = '\0';
																								k=0; t=0; type=0; while (k < j)
																								{
																																while (buffer[k] !=' ' && k < j && buffer[k] !='\0') { block[t] = buffer[k];t++; k++; } k++; 
																																while (buffer[k] ==' ' && k < j) k++; block[t] = '\0'; 
																																if (type <3 ) type++; 
																																else if (type ==3)
																																{
																																								strcpy(Flist[i].list[biter].qv,block);
																																								for (t1=0;t1<Flist[i].list[biter].len;t1++)
																																								{
																																																q = (int)(Flist[i].list[biter].qv[t1])-48; q /=10; q *= -1; 
																																																Flist[i].list[biter].pv[t1] = pow(10,q);
																																																//	fprintf(stdout,"\n PHRED %f \n",pow(10,q)); 
																																								} //fprintf(stderr,"\n");
																																								type =2; biter++;
																																}
																																t=0;
																								} 
																								//fprintf(stdout,"\n%d %s ",Flist[i].blocks,Flist[i].id); 
																								for (j=0;j<Flist[i].blocks;j++) 
																								{
																																//		fprintf(stdout,"| %d %s %s ",Flist[i].list[j].offset,Flist[i].list[j].hap,Flist[i].list[j].qv); 	
																																//for(t1=0;t1<Flist[i].list[j].len;t1++) fprintf(stdout,"%f ",Flist[i].list[j].pv[t1]);
																								}
																								//getchar();
																} fclose(ff);
								}
								/****************************** READ FRAGMENT QUALITY FILE*************************************************/

								struct SNPfrags* snpfrag = (struct SNPfrags*)malloc(sizeof(struct SNPfrags)*snps); for (i=0;i<snps;i++) snpfrag[i].frags = 0;
								// find the first fragment whose endpoint lies at snp 'i' or beyond
								for (i=0;i<snps;i++) snpfrag[i].ff = -1;
								for (i=0;i<fragments;i++)  
								{
																j = Flist[i].list[0].offset; k = Flist[i].list[Flist[i].blocks-1].len + Flist[i].list[Flist[i].blocks-1].offset; 	
																for (t=j;t<k;t++) { if (snpfrag[t].ff == -1) snpfrag[t].ff = i;  } 
								} //for (i=0;i<snps;i++) { fprintf(stdout,"SNP %d firstfrag %d start snp %d \n",i,snpfrag[i].ff,i); } 
								for (i=0;i<fragments;i++)
								{
																Flist[i].cons =1; Flist[i].bad=1;
																for (j=0;j<Flist[i].blocks;j++) 
																{
																								for (k=0;k<Flist[i].list[j].len;k++) snpfrag[Flist[i].list[j].offset+k].frags++; 
																}
								}
								for (i=0;i<snps;i++) { snpfrag[i].flist = (int*)malloc(sizeof(int)*snpfrag[i].frags); snpfrag[i].alist = (char*)malloc(snpfrag[i].frags);}
								update_snpfrags(Flist,fragments,snpfrag,snps,&components);
								for (i=0;i<snps;i++) snpfrag[i].elist = (struct edge*)malloc(sizeof(struct edge)*snpfrag[i].edges);  
								for (i=0;i<snps;i++) snpfrag[i].telist = (struct edge*)malloc(sizeof(struct edge)*snpfrag[i].edges);  
								add_edges(Flist,fragments,snpfrag,snps,&components);
								fprintf(stderr,"fragments %d snps %d component(blocks) %d\n",fragments,snps,components); 

								char* h1 = (char*)malloc(snps+1);		 char* tree_hap = (char*)malloc(snps+1);	
								int* bn = (int*)malloc(4*snps); char* aaron = (char*)malloc(snps+1);
								int* cn = (int*)malloc(4*snps);
								char* tree_besthap = (char*)malloc(snps+1);  char* besthap_mec = (char*)malloc(snps+1);
								for (i=0;i<snps;i++) { h1[i] = '-'; tree_hap[i]='-'; bn[i] = -1; aaron[i] = '-'; cn[i] = -1;} // bn is component no, hold is old haplotype 
								double p=0.02,P=p;				for (i=0;i<snps;i++) snpfrag[i].pv = p;
								time_t ts; time(&ts); srand48((long int)ts);

								/****************************** READ HAPLOTYPE SOLUTION*************************************************/
								FILE* sf = fopen(sol,"r");
								j=0;while (1)
								{
																fscanf(sf,"%s ",id); if (strcmp(id,"BLOCK:") !=0) break; fscanf(sf,"%s %d %s %d %s %d \n",id,&offset,id,&len,id,&phased);  j++;
																for (i=0;i<len;i++) fscanf(sf,"%s %c %c \n",id,&c1,&c2); fscanf(sf,"%s \n",id); 
								}fclose(sf); 

								blocks = j;
								struct BLOCK* blist = (struct BLOCK*)malloc(sizeof(struct BLOCK)*blocks);
								sf = fopen(sol,"r");
								j=0;while (1)
								{
																fscanf(sf,"%s ",id); if (strcmp(id,"BLOCK:") !=0) break;  //fprintf(stdout,"%s %d\n",id,j-1);
																fscanf(sf,"%s %d %s %d %s %d \n",id,&offset,id,&len,id,&phased); blist[j].offset = offset-1; blist[j].length = len; blist[j].phased = phased;
																if (pflag) fprintf(stdout,"BLOCK--- %9d len %5d phased %5d \n",offset,len,phased); j++;
																for (i=0;i<len;i++)
																{
																								fscanf(sf,"%s %c %c \n",id,&c1,&c2); if (c1 != '-') { h1[offset+i-1] = c1; bn[offset+i-1] = offset; aaron[offset+i-1] =c1; snpfrag[offset+i-1].blockno = j;}
																								strcpy(snpfrag[offset+i-1].id,id);
																								// offset is the id of each block since it is supposed to be unique  
																} fscanf(sf,"%s \n",id); 
								}fclose(sf);
								/****************************** READ HAPLOTYPE SOLUTION*************************************************/

								fprintf(stdout,"comparing haplotype file to quality values..... \n");
								//if (QV == 1) compare_Flist_hap(snpfrag,snps,Flist,fragments,aaron,0,QV); //exit(0);

								struct BLOCK* clist = (struct BLOCK*)malloc(sizeof(struct BLOCK)*components); component =0;
								for (i=0;i<snps;i++)	
								{
																if (snpfrag[i].component !=i || snpfrag[i].csize <= 1) continue;
																j=snpfrag[i].csize; t=i;
																while (j > 0 && t < snps) {	if (snpfrag[t].component == snpfrag[i].component) j--; t++; } 
																clist[component].length = t-i; clist[component].phased = snpfrag[i].csize; clist[component].offset = i; 
																clist[component].haplotype = (char*)malloc(t-i+1);  
																for (j=i;j<t;j++) { if (snpfrag[j].component == snpfrag[i].component) cn[j] = i; }
																for (j=i;j<t;j++) { if (snpfrag[j].component == snpfrag[i].component) snpfrag[j].bcomp = component; }
																for (j=i;j<t;j++) { if (snpfrag[j].component == snpfrag[i].component) clist[component].haplotype[j-i] = '0'; else clist[component].haplotype[j-i] = '-';}
																//								fprintf(stdout,"component %d length %d phased %d %d...%d\n",component,clist[component].length,clist[component].phased,clist[component].offset,clist[component].offset+clist[component].length-1);
																//						for (j=0;j<clist[component].length;j++) fprintf(stdout,"%c",clist[component].haplotype[j]); fprintf(stdout,"\n");
																component++;
								} 
								for (i=0;i<components;i++) clist[i].frags=0;			for (i=0;i<fragments;i++) clist[snpfrag[Flist[i].list[0].offset].bcomp].frags++;
								for (i=0;i<components;i++) clist[i].flist = (int*)malloc(4*clist[i].frags); 	for (i=0;i<components;i++) clist[i].frags=0;
								for (i=0;i<fragments;i++) { clist[snpfrag[Flist[i].list[0].offset].bcomp].flist[clist[snpfrag[Flist[i].list[0].offset].bcomp].frags] = i; clist[snpfrag[Flist[i].list[0].offset].bcomp].frags++; } 
								for (i=0;i<components;i++) { clist[i].treecompute ='1'; clist[i].dealloc = '0'; } 

								int tree_bestscore_mec = 0, bestscore_mec = 0,calls=0, miscalls=0,tree_miscalls=0;
								float treebest_ll=0, tree_ll=0, mcmc_ll=0,mcmcbest_ll=0,ll=0; 
								int delta=0,orig=0,block_errors=0, errors=0,Z=0;


								for (i=0;i<components;i++) 
								{ 
																//													if (clist[i].phased > 90 && clist[i].phased < 100) {print_block_frags(clist,i,aaron,h1, cn,Flist,fragments,snpfrag,frags); exit(0); }
								}
								if (SINGLE ==1) {TREE=0; MINCUT =0;}
								if (fraction > 0) ERROR = 1; else ERROR =0; fprintf(stdout,"fraction %d \n",fraction); //getchar();
								if (ERROR ==1)
								{
																for(i=0;i<fragments;i++) correct_fragment(Flist,i,aaron);  errors= mutate_Flist(Flist,fragments,(double)fraction/(double)1000);
																fprintf(stderr,"# base calls flipped %d \n",errors);
																fprintf(stdout,"# base calls flipped %d \n",errors);  getchar();
																AARON =0; QV = -1;
																//print_block_frags(clist,0,aaron,aaron, cn,Flist,fragments,snpfrag,frags); exit(0); 
																TREE = 0; 
								}

								mecscore(Flist,fragments,aaron,&ll,&calls,&miscalls);
								fprintf(stdout,"input haplotype MEC %d calls %d log likelihood %f\n",miscalls,calls,ll); //getchar();
        // for ERROR =1, this is just measuring the edit distance of mutated matrix from perfect haplotype which is expected
        // to be of the same order as the # of flips performed 
								float errprob = (double)miscalls/(double)calls;  if (EXAMPLE) errprob = 0.05;

								if (STDERR) fprintf(stderr,"input haplotype MEC %d calls %d log likelihood %f\n",miscalls,calls,ll); 
																if (EXAMPLE) getchar();

								if (QV != 1)	
								{
																for (i=0;i<fragments;i++) {for (j=0;j<Flist[i].blocks;j++){  for (k=0;k<Flist[i].list[j].len;k++) Flist[i].list[j].pv[k] = errprob;  } }
																fprintf(stderr,"global q value %f \n",(double)miscalls/(double)calls); //exit(0);
								}

								if (AARON ==0)        
								{	
																for (i=0;i<snps;i++) h1[i] = '-';
																fprintf(stdout,"USING fragment clustering to obtain an INITIAL SOLUTION \n");
																frag_cluster_initialize(Flist,fragments,snpfrag,h1,snps,clist,components);
																if (RANDOM_START)
																{				
																								for (i=0;i<snps;i++) { if (h1[i] != '-') { if (drand48() < 0.5) h1[i] = '0'; else  h1[i] = '1'; } } 
																								//								for (i=0;i<snps;i++) { if (h1[i] != '-')  h1[i] = '0';} 
																}
																for (i=0;i<snps;i++) tree_hap[i] = h1[i]; 
																for (i=0;i<snps;i++) { besthap_mec[i] = h1[i]; tree_besthap[i] =h1[i];} 
								}
								else 
								{
																for (i=0;i<snps;i++) h1[i] = aaron[i]; 
																for (i=0;i<snps;i++) tree_hap[i] = aaron[i]; 
																for (i=0;i<snps;i++) { besthap_mec[i] = aaron[i]; tree_besthap[i] =aaron[i];} 
								}

								// for each block, we maintain best haplotype solution under MFR criterion 
								// compute the component-wise score for 'aaron' haplotype 
								miscalls=0;bestscore_mec=0; mcmcbest_ll= mcmc_ll = 0;
								for (k=0;k<components;k++)
								{
																clist[k].MEC =0; clist[k].bestMEC =0; clist[k].calls =0;	clist[k].treeMEC =0; clist[k].treebestMEC =0; 
																clist[k].LL = 0; clist[k].treeLL =0;
																for (i=0;i<clist[k].frags;i++) 
																{
																								update_fragscore(Flist,clist[k].flist[i],h1); 
																								clist[k].MEC += Flist[clist[k].flist[i]].currscore;
																								clist[k].LL += Flist[clist[k].flist[i]].ll;
																								clist[k].calls += Flist[clist[k].flist[i]].calls;
																								clist[k].treeMEC += compute_fragscore(Flist,clist[k].flist[i],tree_besthap,&ll);  clist[k].treeLL += ll; 

																} 
																clist[k].bestMEC = clist[k].MEC; bestscore_mec += clist[k].bestMEC; miscalls += clist[k].MEC;	
																clist[k].bestLL = clist[k].LL; mcmcbest_ll += clist[k].bestLL; mcmc_ll += clist[k].LL;	
																clist[k].treebestLL = clist[k].treeLL; 
																clist[k].lastMEC = clist[k].bestMEC; clist[k].treebestMEC = clist[k].treeMEC; 
																fprintf(stdout,"component %d offset %d length %d phased %d  calls %d MEC %d LL %f\n",k,clist[k].offset,clist[k].length,clist[k].phased,clist[k].calls,clist[k].MEC,clist[k].LL);

								}


								for (i=0;i<snps;i++) snpfrag[i].island = '0';
								for (i=0;i<fragments;i++) {for (j=0;j<Flist[i].blocks;j++){  for (k=0;k<Flist[i].list[j].len-1;k++) { snpfrag[Flist[i].list[j].offset+k].island = '1'; } } }

								slist = (int*)malloc(sizeof(int)*snps);
								float rd1 = drand48();
								int maxiter=500000; if (EXAMPLE ==1) maxiter = 100000000;
								char mincut_temp[200]; sprintf(mincut_temp,"%s.mincut.output-%f",frags,rd1);
								for (iter=0;iter<maxiter;iter++)
								{
																for (k=0;k<components;k++) // COMPUTATION OF TREE FOR EACH COMPONENT 
																{
																								//if (clist[k].phased > 100) continue;
																								if (clist[k].treecompute =='1') 
																								{
																																i=0;for (j=clist[k].offset;j<clist[k].offset+clist[k].length;j++) 
																																{
																																								if (snpfrag[clist[k].offset].component == snpfrag[j].component) { slist[i] = j; i++; } 
																																}
																																clist[k].nodes=0;
																																fprintf(stdout,"\n component %d length %d phased %d %d...%d \n",k,clist[k].length,clist[k].phased,clist[k].offset,clist[k].offset+clist[k].length-1);
																																if (STDERR) fprintf(stderr,"\n component %d length %d phased %d %d...%d \n",k,clist[k].length,clist[k].phased,clist[k].offset,clist[k].offset+clist[k].length-1);
																																if (MINCUT ==1 && SINGLE==0) { MCout = fopen(mincut_temp,"w"); partitions(snpfrag,snps,slist,clist[k].phased,min,Flist,fragments,MCout,&clist[k].nodes,h1); fclose(MCout);}
																																else { MCout = fopen(mincut_temp,"w"); NJtree(snpfrag,snps,clist,k,MCout,&clist[k].nodes); fclose(MCout); } 
																																fprintf(stdout,"\n completed MIN-CUT for component %d length %d phased %d %d...%d tree nodes %d\n",k,clist[k].length,clist[k].phased,clist[k].offset,clist[k].offset+clist[k].length-1,clist[k].nodes);
																																if (STDERR) fprintf(stderr,"\n completed MIN-CUT for component %d length %d phased %d %d...%d tree nodes %d\n",k,clist[k].length,clist[k].phased,clist[k].offset,clist[k].offset+clist[k].length-1,clist[k].nodes);
																																clist[k].tree = (struct NODE*)malloc(sizeof(struct NODE)*clist[k].nodes);
																																MCout = fopen(mincut_temp,"r"); 
																																for (i=0;i<clist[k].nodes;i++)
																																{
																																								fscanf(MCout,"%d %d %d ",&clist[k].tree[i].offset,&clist[k].tree[i].length,&clist[k].tree[i].phased); 
																																								clist[k].tree[i].flips =0; clist[k].tree[i].potflips =0; clist[k].tree[i].updatenodes=-1;
																																								clist[k].tree[i].p00 = 100; clist[k].tree[i].p01 = 100;

																																								if (clist[k].tree[i].phased <= min) clist[k].tree[i].leaf = '1'; else clist[k].tree[i].leaf = '0'; // leaf 
																																								clist[k].tree[i].parent = clist[k].tree[i].ch1 = clist[k].tree[i].ch2 = -1;
																																								clist[k].tree[i].bitvec = (char*)malloc(clist[k].tree[i].length+1);
																																								for (j=0;j<clist[k].tree[i].length;j++) fscanf(MCout,"%c ",&clist[k].tree[i].bitvec[j]); 
																																								fscanf(MCout,"\n");

																																								clist[k].tree[i].frags = -1; compute_flist(snpfrag,snps,Flist,fragments,clist,k,i,0);
																																								if (clist[k].tree[i].frags > 0) clist[k].tree[i].flist = (int*)malloc(4*clist[k].tree[i].frags);
																																								compute_flist(snpfrag,snps,Flist,fragments,clist,k,i,0);

																																								clist[k].tree[i].sfrags = -1; compute_flist(snpfrag,snps,Flist,fragments,clist,k,i,1);
																																								if (clist[k].tree[i].sfrags > 0) clist[k].tree[i].sflist = (int*)malloc(4*clist[k].tree[i].sfrags);
																																								compute_flist(snpfrag,snps,Flist,fragments,clist,k,i,1);
																																}
																																fclose(MCout);
																																for (i=0;i<clist[k].nodes;i++)
																																{
																																								if (UPDATE ==0) continue; 
																																								compute_updatelist(snpfrag,snps,Flist,fragments,clist,k,i); 
																																								if (clist[k].tree[i].updatenodes >0) clist[k].tree[i].updatelist = (int*)malloc(4*clist[k].tree[i].updatenodes); 
																																								compute_updatelist(snpfrag,snps,Flist,fragments,clist,k,i);

																																} 

																																if (MINCUT ==1) construct_tree(clist,k,0,1); 
																																clist[k].treecompute ='0'; 
																								}
																								for (t1=0;t1<clist[k].nodes;t1++)
																								{
																																t = (int)((drand48()-0.0000001)*clist[k].nodes); flip=0;	
																																if (clist[k].tree[t].phased >1 && SINGLE ==1) continue;
																																if ((clist[k].tree[t].phased != 1 && clist[k].tree[t].phased != snps/2) && EXAMPLE ==1) continue;
																																flip = sample_hap_rest(snpfrag,snps,Flist,fragments,clist,k,t,h1); blockflips++;  
																																if (flip ==1) 
																																{
																																								if (clist[k].tree[t].phased <= min && min > 1) realallflips++; else realblockflips++; 
																																								bestscore_mec -= clist[k].bestMEC; miscalls -= clist[k].MEC;
																																								mcmcbest_ll -= clist[k].bestLL;  mcmc_ll -= clist[k].LL;	
																																								for (i=0;i<clist[k].tree[t].frags;i++)
																																								{
																																																j = clist[k].tree[t].flist[i];
																																																clist[k].LL -= Flist[j].ll; 	 clist[k].MEC -= Flist[j].currscore;
																																																update_fragscore(Flist,j,h1);							
																																																clist[k].LL += Flist[j].ll; clist[k].MEC += Flist[j].currscore;

																																								}
																																								//if (clist[k].MEC < clist[k].bestMEC) clist[k].bestMEC = clist[k].MEC;
																																								if (clist[k].LL > clist[k].bestLL+0.01)
																																								{
																																																for (i=0;i<clist[k].length;i++) 
																																																{
																																																								if (clist[k].haplotype[i] !='-') besthap_mec[clist[k].offset+i] = h1[clist[k].offset+i];
																																																} clist[k].bestLL = clist[k].LL;
																																																clist[k].bestMEC = clist[k].MEC;
																																								}
																																								bestscore_mec += clist[k].bestMEC; 		miscalls+= clist[k].MEC;
																																								mcmcbest_ll += clist[k].bestLL;  mcmc_ll += clist[k].LL;
																																}
																								} 
																								if ( (iter%2000==0 && iter > 0 && iter < 5000 ))// || iter > 10000 && iter%10000==0)
																								{ 
																																if (SINGLE ==0 && UPDATETREE ==1) { clist[k].treecompute ='1'; clist[k].dealloc ='1'; clist[k].lastMEC = clist[k].bestMEC;	}
																								} 
																								if (clist[k].dealloc =='1') 
																								{
																																for (i=0;i<clist[k].nodes;i++) free(clist[k].tree[i].bitvec);
																																for (i=0;i<clist[k].nodes;i++) { if (clist[k].tree[i].frags >0) free(clist[k].tree[i].flist); }
																																for (i=0;i<clist[k].nodes;i++) { if (clist[k].tree[i].sfrags >0) free(clist[k].tree[i].sflist); }
																																if (UPDATE ==1)
																																{
																																								for (i=0;i<clist[k].nodes;i++) { if (clist[k].tree[i].updatenodes >0) free(clist[k].tree[i].updatelist); } 
																																}
																																free(clist[k].tree);  clist[k].dealloc = '0'; 
																								}
																								if ((iter%1000 ==0 && GIBBS ==1 && EXAMPLE ==0) || (iter%(snps*snps*((int)log2(snps))) ==0 && EXAMPLE ==1 && iter > 100000))
																								{
//																																if (STDERR) fprintf(stderr,"component %d offset %d phased %d length %d MCMC %6d %d %d %f BESTMEC %d Log-lh %f bestLL %f\n",k,clist[k].offset,clist[k].phased,clist[k].length,iter,clist[k].MEC,clist[k].calls,(double)clist[k].MEC/(double)clist[k].calls,clist[k].bestMEC,clist[k].LL,clist[k].bestLL);
																																fprintf(stdout,"component %d offset %d phased %d length %d MCMC %6d %d %d %f BESTMEC %d Log-lh %f bestLL %f\n",k,clist[k].offset,clist[k].phased,clist[k].length,iter,clist[k].MEC,clist[k].calls,(double)clist[k].MEC/(double)clist[k].calls,clist[k].bestMEC,clist[k].LL,clist[k].bestLL);
//																																if (clist[k].MEC == MF && EXAMPLE ==1) { fprintf(stdout,"sampled best haplotype iter %d \n",iter);} // exit(0);}
																																if (EXAMPLE ==1) { for (i=0;i<snps;i++) fprintf(stdout,"%c",h1[i]); if (h1[snps/2] == h1[snps/2-1] && clist[k].MEC == MF ) fprintf(stdout," H1\n"); else if (clist[k].MEC == MF) fprintf(stdout," H2\n"); else fprintf(stdout,"\n");} 
	//																															if (EXAMPLE ==1) { for (i=0;i<snps;i++) fprintf(stderr,"%c",h1[i]); if (h1[snps/2] == h1[snps/2-1] && clist[k].MEC == MF ) fprintf(stderr," H1\n"); else if (clist[k].MEC == MF) fprintf(stderr," H2\n"); else fprintf(stderr,"\n");} 
																																//if (clist[k].bestMEC == 0 && iter >= 1000) maxiter = 1000;

																								}
																}
																if ((iter%10 ==0 && GIBBS ==1 && EXAMPLE ==0))
																{
																								if (STDERR) fprintf(stderr,"MCMC %6d %d %d %f BESTMEC %d Log-lh %f bestLL %f\n",iter,miscalls,calls,(double)miscalls/(double)calls,bestscore_mec,mcmc_ll,mcmcbest_ll);
																								fprintf(stdout,"FULLMCMC %6d %d %d %f BESTMEC %d Log-lh %f bestLL %f\n",iter,miscalls,calls,(double)miscalls/(double)calls,bestscore_mec,mcmc_ll,mcmcbest_ll);
																								if (iter%20000==0 && iter > 0 && ERROR ==0 && SINGLE ==0) { fprintf(stdout,"iter %d OUTPUTTING BEST SOLUTION BY MEC USING ALL FRAGMENTS \n",iter); print_blocks(clist,components,tree_besthap,besthap_mec,h1,cn,Flist,fragments,bestscore_mec,snpfrag);}
																								if (iter%100000==0 && iter > 0 && ERROR ==0 && SINGLE ==1) { fprintf(stdout,"iter %d OUTPUTTING BEST SOLUTION BY MEC USING ALL FRAGMENTS \n",iter); print_blocks(clist,components,tree_besthap,besthap_mec,h1,cn,Flist,fragments,bestscore_mec,snpfrag);}
																								if (iter%1000==0 && iter > 0 && ERROR ==1)
																								{
																																block_errors = compare_haps(clist,components,aaron,besthap_mec,snpfrag,snps);
																																fprintf(stderr,"errors in phasing %d \n",block_errors);
																																fprintf(stdout,"errors in phasing %d \n",block_errors);
																								}

																}
																if (iter >= burnin && iter%thinrate ==0 && EXAMPLE ==0)
																{
																								if (Z%200 ==0 && Z > 0 )
																								{
																																fprintf(stderr,"OUTPUTTING posterior error probabilities for base calls \n");
																																fprintf(stdout,"OUTPUTTING posterior error probabilities for base calls \n");
																																compare_Flist_hap(snpfrag,snps,Flist,fragments,aaron,Z,QV);
																								}
																								Z++; update_posterior(Flist,fragments,h1,Z);
																}
								}
}


