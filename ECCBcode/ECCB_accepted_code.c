// CODE STARTED SEPT 10 2007 4pm //  

// april 8 2008 this code used for producing results in ECCB 2008 paper

#include "functions_ECCB.h"
int TREE =0; 
int AARON =0;
int STDERR =1;
int RANDOM_START=1;
int USE_HAP = 1;
int ERROR = 1;
int QV = -1;
int MINCUT =1;
int burnin = 20000;
int thinrate = 100;
int MCMCruns=0;
int SINGLE =0;
int UPDATE =1;
int MINCUTALGO =1;


int fraglength(struct fragment* Flist,int f)
{
								int i=0,l=0; 
								for (i=0;i<Flist[f].blocks;i++) l += Flist[f].list[i].len;  
								return l; 
}

int main(int argc, char** argv)
{
								if (argc < 2) { fprintf(stderr," enter some arguments ..... exiting \n"); exit(0);}
								if (MINCUTALGO ==2) RANDOM_START=0;

								char frags[100]; char sol[100]; 
								int flip =0;			if (argc >= 4) flip = atoi(argv[3]);
								if (argc <3) return 1;
								maxcut_opt(argv[1],argv[2],flip); exit(0);
}

/**************** DETERMINISTIC MAX_CUT MEC IMPLEMENTATION *********************************************************//////


float edge_weight(char* hap,int i, int j, char* p)
{	
								if (hap[i] == hap[j] && p[0] == p[1]) return 1;
								if (hap[i] != hap[j] && p[0] != p[1]) return 1;
								if (hap[i] == hap[j] && p[0] != p[1]) return -1;
								if (hap[i] != hap[j] && p[0] == p[1]) return -1;
}

float compute_goodcut(struct SNPfrags* snpfrag,char* hap,int* slist,int N,struct fragment* Flist, int algo,int iteration)
{
								// given a haplotype 'hap' and a fragment matrix, find a cut with positive score 
								int totaledges=0,i=0,j=0,k=0,l=0,fmec=0,pflag=0;  int wf = 0; //if (drand48() < 0.5) wf=1;
        float W=0, Wpos =0,Wneg=0,a,b,ll;
								for (i=0;i<N;i++)
								{
																snpfrag[slist[i]].tedges=0; k=-1;
																for (j=0;j<snpfrag[slist[i]].edges;j++) 
																{
																								if (k != snpfrag[slist[i]].elist[j].snp) { snpfrag[slist[i]].tedges++; k = snpfrag[slist[i]].elist[j].snp; } 
																}
								}
								for (i=0;i<N;i++)
								{
																snpfrag[slist[i]].tedges=0; k=-1;
																for (j=0;j<snpfrag[slist[i]].edges;j++) 
																{
																								if (k != snpfrag[slist[i]].elist[j].snp) 
																								{ 
																																snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].snp = snpfrag[slist[i]].elist[j].snp; 
																																k = snpfrag[slist[i]].elist[j].snp; 
																																W = (float)edge_weight(hap,slist[i],k,snpfrag[slist[i]].elist[j].p);
																																if (wf ==0) W /= (fraglength(Flist,snpfrag[slist[i]].elist[j].frag)-1);	
																																snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].w = W; 
																																snpfrag[slist[i]].tedges++;  totaledges++;
																								} 
																								else if (k == snpfrag[slist[i]].elist[j].snp) 
																								{
																																W = (float)edge_weight(hap,slist[i],k,snpfrag[slist[i]].elist[j].p);
																																if (wf ==0) W /= (fraglength(Flist,snpfrag[slist[i]].elist[j].frag)-1); 
																																snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges-1].w += W; 
																								}
																}
								}
								// edge contraction algorithm: merge vertices until only two nodes left or total edge weight of graph is negative  
								int startnode =(int)(drand48()*N); if (startnode ==N) startnode--;  
								int secondnode= -1 ;  // root of 2nd cluster initially not there
								// chose a positive edge to initialize the two clusters and run this algorithm $O(m)$ times for each block 
								// a negative weight cut should have at least one negative edge or if there is no negative weight edge, the edge with lowest weight 
								int* score1= (int*)malloc(sizeof(int)*N);  
								int* score2= (int*)malloc(sizeof(int)*N);  

								for (i=0;i<N;i++) snpfrag[slist[i]].parent = slist[i];  int V = N, Asize; 
								float bscore,curr_cut=0,best_cut=10000;
								int last1,last2,best_point,snp_add;
								int moved =1,c1=0,c2=0; char* mincut; char* bestmincut;
								int size_small,best_small=0,secondlast=0,last=0,iter=0,maxiter=N/10; if (N/10 < 1) maxiter =1; //maxiter = 1;

								/*****************************Maintain two clusters and add each vertex to one of these two ******************/
								if (algo ==1)
								{
																bestmincut = (char*)malloc(N); for (i=0;i<N;i++) bestmincut[i] = '0';
																//for (iter=0;iter<totaledges*(int)(log2(totaledges));iter++)
																for (iter=0;iter<maxiter;iter++)
																{
																								i = (int)(drand48()*totaledges-0.0001); j=0;
																								while (i >= snpfrag[slist[j]].tedges) { i -= snpfrag[slist[j]].tedges; j++;} 
																								startnode = slist[j]; secondnode = snpfrag[slist[j]].telist[i].snp; 
																								if (snpfrag[slist[j]].telist[i].w >=1) continue; 
																								for (i=0;i<N;i++) snpfrag[slist[i]].parent = slist[i];
																								V = N; 
																								while (V > 2) // more than two clusters 
																								{
																																bscore =-1000; snp_add = -1; 
																																for (i=0;i<N;i++) 
																																{
																																								snpfrag[slist[i]].score  = 0;
																																								if (snpfrag[slist[i]].parent == startnode) continue;
																																								if ( snpfrag[slist[i]].parent == secondnode) continue; 
																																								for (j=0;j<snpfrag[slist[i]].tedges;j++) 
																																								{
																																																if (snpfrag[snpfrag[slist[i]].telist[j].snp].parent == startnode) {snpfrag[slist[i]].score += snpfrag[slist[i]].telist[j].w; score1[i] += snpfrag[slist[i]].telist[j].w;} 
																																																if (snpfrag[snpfrag[slist[i]].telist[j].snp].parent == secondnode) {snpfrag[slist[i]].score -= snpfrag[slist[i]].telist[j].w; score2[i] += snpfrag[slist[i]].telist[j].w;}
																																								}
																																								if (snpfrag[slist[i]].score > bscore) { bscore = snpfrag[slist[i]].score; snp_add = i; } 
																																								else if (-1*snpfrag[slist[i]].score > bscore) { bscore = -1*snpfrag[slist[i]].score; snp_add = i; } 
																																}
																																if (bscore != -1000)
																																{
																																								if (snpfrag[slist[snp_add]].score > 0) { snpfrag[slist[snp_add]].parent = startnode; V--; } 
																																								else if (snpfrag[slist[snp_add]].score < 0 ) {snpfrag[slist[snp_add]].parent = secondnode; V--; } 
																																								else if (snpfrag[slist[snp_add]].score ==0) 
																																								{
																																																if (drand48() < 0.5) snpfrag[slist[snp_add]].parent = startnode; else snpfrag[slist[snp_add]].parent = secondnode; V--; 
																																								}
																																}

																								}
																								// compute score of the cut computed above 
																								for (i=0;i<N;i++) { if (snpfrag[slist[i]].parent == startnode) snpfrag[slist[i]].parent =0; else snpfrag[slist[i]].parent = 1;  } 
																								c1=0;c2=0; for (i=0;i<N;i++) { if (snpfrag[slist[i]].parent == 0) c1++; else c2++;  } 
																								curr_cut =0; 
																								for(i=0;i<N;i++) 
																								{ 
																																for (j=0;j<snpfrag[slist[i]].tedges;j++)
																																{
																																								if (snpfrag[snpfrag[slist[i]].telist[j].snp].parent != snpfrag[slist[i]].parent) 		curr_cut += snpfrag[slist[i]].telist[j].w/2; 
																																}
																								}
																								moved =1; while (moved > 0) // any improvement in score of cut  
																								{
																																moved =0;
																																for (i=0;i<N;i++)
																																{
																																								snpfrag[slist[i]].score  = 0;
																																								for (j=0;j<snpfrag[slist[i]].tedges;j++)
																																								{
																																																if (snpfrag[snpfrag[slist[i]].telist[j].snp].parent == 0) snpfrag[slist[i]].score += snpfrag[slist[i]].telist[j].w;
																																																if (snpfrag[snpfrag[slist[i]].telist[j].snp].parent == 1) snpfrag[slist[i]].score -= snpfrag[slist[i]].telist[j].w;
																																								}
																																								if (snpfrag[slist[i]].parent ==0 && snpfrag[slist[i]].score < 0 && c1 > 1)
																																								{
																																																snpfrag[slist[i]].parent =1; curr_cut += snpfrag[slist[i]].score; moved++; c1--; c2++;
																																								}
																																								if (snpfrag[slist[i]].parent ==1 && snpfrag[slist[i]].score > 0 && c2 > 1)
																																								{
																																																snpfrag[slist[i]].parent =0; curr_cut -= snpfrag[slist[i]].score; moved++; c2--; c1++;
																																								}
																																}
																								} if (c1 ==0 || c2==0) { fprintf(stdout," cut size is 0 red \n"); exit(0); } 
																								if (curr_cut < best_cut)
																								{
																																best_cut = curr_cut; 
                  for (i=0;i<N;i++) { if (snpfrag[slist[i]].parent == 1) bestmincut[i] ='1'; else bestmincut[i] = '0'; } 
																								}
//																								if (best_cut < -5) iter = maxiter; 
																}
																for(i=0;i<N;i++) { if (bestmincut[i] == '1') slist[i]= -1*slist[i]-1; }
																free(score1); free(score2); free(bestmincut);
																return best_cut;
								}
								/*
											fprintf(stdout,"startnode %d secondnode %d %f \n",slist[startnode],slist[secondnode],best_cut);
											for (i=0;i<N;i++)
											{
											fprintf(stdout,"\ncol %d edges %d par %d| ",slist[i],snpfrag[slist[i]].tedges,snpfrag[slist[i]].parent);
											for (j=0;j<snpfrag[slist[i]].tedges;j++) 
											{
											fprintf(stdout,"(%d %f) ",snpfrag[slist[i]].telist[j].snp,snpfrag[slist[i]].telist[j].w); 
											}
											} getchar();
											*/

								/*****************************Maintain two clusters and add each vertex to one of these two ******************/

								/************************ NOW apply SIMPLE MIN CUT algorithm of wagner ESA 1994 ********************************/
								if (algo ==2)
								{
																mincut = (char*)malloc(N); bestmincut = (char*)malloc(N); 
																for (i=0;i<N;i++) bestmincut[i] = '0';
																while (V >1)
																{
																								Asize = 1;
																								for (i=0;i<N;i++) snpfrag[slist[i]].Aset = '0';  snpfrag[slist[0]].Aset = '1'; 
																								for (i=0;i<N;i++) mincut[i] = '0'; mincut[0] = '1';
																								for (i=0;i<N;i++) snpfrag[slist[i]].score  = 0; bscore =-1000;
																								for (j=0;j<snpfrag[slist[startnode]].tedges;j++)
																								{
																																if (snpfrag[snpfrag[slist[startnode]].telist[j].snp].Aset =='0') snpfrag[snpfrag[snpfrag[slist[startnode]].telist[j].snp].parent].score += snpfrag[slist[startnode]].telist[j].w;
																								}
																								while (Asize < V)
																								{
																																bscore = -1000;
																																for (i=0;i<N;i++) { if (snpfrag[slist[i]].score > bscore && snpfrag[slist[i]].Aset == '0') {bscore = snpfrag[slist[i]].score; snp_add = i; } }
																																for (i=0;i<N;i++) { if (snpfrag[slist[i]].parent == slist[snp_add]) snpfrag[slist[i]].Aset = '1'; }

																																for (i=0;i<N;i++) 
																																{
																																								if (snpfrag[slist[i]].parent == slist[snp_add])
																																								{
																																																for (j=0;j<snpfrag[slist[i]].tedges;j++)
																																																{
																																																								if (snpfrag[snpfrag[slist[i]].telist[j].snp].Aset =='0') snpfrag[snpfrag[snpfrag[slist[i]].telist[j].snp].parent].score += snpfrag[slist[i]].telist[j].w;
																																																}
																																																if (Asize <V-1) mincut[i] = '1';
																																								}
																																}
																																Asize++;
																																if (pflag) fprintf(stdout,"(%d %f)",slist[snp_add],bscore);
																																if (Asize == V-1) secondlast = snp_add;  if (Asize == V) last = snp_add;
																								}
																								// also maintain last two vertices added and they should be merged for next round of mincut computation 
																								if (slist[secondlast] > slist[last]) snpfrag[slist[secondlast]].parent = snpfrag[slist[last]].parent; 
																								else  snpfrag[slist[last]].parent = snpfrag[slist[secondlast]].parent;
																								for (i=0;i<N;i++) snpfrag[slist[i]].parent =  snpfrag[snpfrag[slist[i]].parent].parent;
																								curr_cut = bscore;  
																								size_small =0; for (i=0;i<N;i++) {if (mincut[i] == '0') size_small++; } if (size_small >= N/2) size_small = N-size_small;
																								for (i=0;i<N;i++) snpfrag[slist[i]].parent =  snpfrag[snpfrag[slist[i]].parent].parent;
																								if (curr_cut < best_cut || (size_small > best_small && curr_cut <= best_cut )) 
																								{ 
																																for (i=0;i<N;i++) bestmincut[i]= mincut[i];  best_small = size_small;   
																																best_cut = curr_cut; best_point = V; 
																								}
																								V--;
																} 
																// compute the actual min cut 
																if (pflag) fprintf(stdout,"MIN-CUT value %d point %d \n",best_cut,best_point); 
																if (pflag) for (i=0;i<N;i++) fprintf(stdout,"%c",bestmincut[i]) ;
																if (pflag) { fprintf(stdout,"BLOCK1 "); for (i=0;i<N;i++) { if (bestmincut[i] == '0') fprintf(stdout,"%d ",slist[i]);} fprintf(stdout,"\n");}
																if (pflag) { fprintf(stdout,"BLOCK2 "); for (i=0;i<N;i++) { if (bestmincut[i] == '1') fprintf(stdout,"%d ",slist[i]);} fprintf(stdout,"\n");}
																for(i=0;i<N;i++) { if (bestmincut[i] == '1') slist[i]= -1*slist[i]-1; }
																free(score1); free(score2);
																free(mincut); free(bestmincut); return best_cut;
								}
								/************************ NOW apply SIMPLE MIN CUT algorithm of wagner ESA 1994 ********************************/
}


int maxcut_opt(char* frags,char* sol,int fraction)
{
								char newfrags[200]; char qvfile[200];
								char command[500]; char command1[500];
								sprintf(command,"more %s | awk 'BEGIN {rows=1; snps=0;} { if ($3 == \"\") snps = $2; else if ($4 != NULL || length($3) > 1) {rows =rows+1; print $0;} } END { print rows,snps;} ' |  awk ' { flag =0; f =4; for (i =3; i<= 100; i+=1) {  if ($i == NULL && flag ==0  && $i != 0) { f= i; flag =1; }}  t = f-2; if (f >=4) print t/2,$0;  else print $0; }' | sort -g -k 3 -k 1 > %s.SORTED ;",frags,frags);
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
								for (i=0;i<snps;i++) { h1[i] = '0'; tree_hap[i]='0'; bn[i] = -1; aaron[i] = '0'; cn[i] = -1;} // bn is component no, hold is old haplotype 
								double p=0.02,P=p;				for (i=0;i<snps;i++) snpfrag[i].pv = p;
								time_t ts; time(&ts); srand48((long int)ts);

								struct tm  *ts1;   char       buf[80];								time_t     now;		
								time(&now); ts1 = localtime(&now); strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", ts1); printf("%s\n", buf); 

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
								float cutvalue =0;

								mecscore(Flist,fragments,aaron,&ll,&calls,&miscalls);
								fprintf(stdout,"input haplotype MEC %d calls %d log likelihood %f avg fragment length %f\n",miscalls,calls,ll,(float)calls/(float)fragments); //getchar();
								// for ERROR =1, this is just measuring the edit distance of mutated matrix from perfect haplotype which is expected
								// to be of the same order as the # of flips performed 
								float errprob = (double)miscalls/(double)calls;//  if (EXAMPLE) errprob = 0.05;

								if (STDERR) fprintf(stderr,"input haplotype MEC %d calls %d log likelihood %f\n",miscalls,calls,ll); 
								if (QV != 1)	
								{
																for (i=0;i<fragments;i++) {for (j=0;j<Flist[i].blocks;j++){  for (k=0;k<Flist[i].list[j].len;k++) Flist[i].list[j].pv[k] = errprob;  } }
																fprintf(stderr,"global q value %f \n",(double)miscalls/(double)calls); //exit(0);
								}

								if (AARON ==0)        
								{	
																if (RANDOM_START ==1)
																{				
																								fprintf(stdout,"starting from a completely random solution SOLUTION \n");
																								for (i=0;i<snps;i++) { if (h1[i] != '-') { if (drand48() < 0.5) h1[i] = '0'; else  h1[i] = '1'; } } 
																}
																else
																{
																								for (i=0;i<snps;i++) h1[i] = '-';
																								fprintf(stdout,"USING fragment clustering to obtain an INITIAL SOLUTION \n");
																								frag_cluster_initialize(Flist,fragments,snpfrag,h1,snps,clist,components);
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
																clist[k].MEC =0; clist[k].bestMEC =0; clist[k].calls =0;	clist[k].treeMEC =0; clist[k].treebestMEC =0; clist[k].lastMEC = 0; 
																clist[k].LL = 0; clist[k].treeLL =0;
																for (i=0;i<clist[k].frags;i++) 
																{
																								update_fragscore(Flist,clist[k].flist[i],h1); 
																								clist[k].MEC += Flist[clist[k].flist[i]].currscore;

																								update_fragscore(Flist,clist[k].flist[i],aaron); 
																								clist[k].lastMEC += Flist[clist[k].flist[i]].currscore;

																								clist[k].LL += Flist[clist[k].flist[i]].ll;
																								clist[k].calls += Flist[clist[k].flist[i]].calls;
																								clist[k].treeMEC += compute_fragscore(Flist,clist[k].flist[i],tree_besthap,&ll);  clist[k].treeLL += ll; 

																} 
																clist[k].bestMEC = clist[k].MEC; 
																bestscore_mec += clist[k].bestMEC; miscalls += clist[k].MEC;	
																clist[k].bestLL = clist[k].LL; mcmcbest_ll += clist[k].bestLL; mcmc_ll += clist[k].LL;	
																clist[k].treebestLL = clist[k].treeLL; 
								}

								slist = (int*)malloc(sizeof(int)*snps); int tit=0,mecdelta=0;
								int maxiter=100;
								for (iter=0;iter<= maxiter;iter++)
								{
																mecscore(Flist,fragments,h1,&ll,&calls,&miscalls);
        								time(&now); ts1 = localtime(&now); strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", ts1); 
																fprintf(stdout,"iter %d input hap MEC %d calls %d log likelihood %f %s \n",iter,miscalls,calls,ll,buf);  
																fprintf(stderr,"iter %d input hap MEC %d calls %d log likelihood %f %s \n",iter,miscalls,calls,ll,buf);  
																if (iter == maxiter) break;
																for (k=0;k<components;k++) // COMPUTATION OF TREE FOR EACH COMPONENT 
																{
                        if (k%100 ==0) fprintf(stderr,"#");
//																								if (k != 433) continue;
																								if (iter ==0)			
																								{
																														fprintf(stdout,"\n component %d length %d phased %d %d...%d \n",k,clist[k].length,clist[k].phased,clist[k].offset,clist[k].offset+clist[k].length-1); 
										//																				for (j=0;j<clist[k].frags;j++)
								//																						{
		//																																fprintf(stdout," \n frag %d ",clist[k].flist[j]);
				//																														for (t=0;t<Flist[clist[k].flist[j]].blocks;t++) fprintf(stdout," %d %s ",Flist[clist[k].flist[j]].list[t].offset,Flist[clist[k].flist[j]].list[t].hap);  
						//																								} getchar();
																								}
																								i=0;for (j=clist[k].offset;j<clist[k].offset+clist[k].length;j++) 
																								{
																																if (snpfrag[clist[k].offset].component == snpfrag[j].component) { slist[i] = j; i++; } 
																								}
																								// run gamma_1 algorithm, if no improvement using gamma1 then use bigger cuts 
																								for (t=0;t<clist[k].phased;t++)
																								{
																																if (h1[slist[t]] == '1') h1[slist[t]] = '0';
																																else if (h1[slist[t]] == '0') h1[slist[t]] = '1';
																																clist[k].bestMEC = clist[k].MEC; clist[k].MEC =0; 
																																for (i=0;i<clist[k].frags;i++) 
																																{
																																								update_fragscore(Flist,clist[k].flist[i],h1); 
																																								clist[k].MEC += Flist[clist[k].flist[i]].currscore;
																																}
																																if (clist[k].MEC > clist[k].bestMEC) 
																																{
																																								if (h1[slist[t]] == '1') h1[slist[t]] = '0';
																																								else if ( h1[slist[t]] == '0') h1[slist[t]] = '1';
																																								clist[k].MEC = clist[k].bestMEC; 
																																}
																								}
																								i=0;for (j=clist[k].offset;j<clist[k].offset+clist[k].length;j++) 
																								{
																																if (snpfrag[clist[k].offset].component == snpfrag[j].component) { slist[i] = j; i++; } 
																								}

																								cutvalue =10;
																								if (clist[k].MEC > 0) cutvalue = compute_goodcut(snpfrag,h1,slist,clist[k].phased,Flist,MINCUTALGO,iter);
																								// flip the subset of columns in slist with positive value 
																								if (cutvalue <= 3 || MINCUTALGO ==2) 
																								{	 //getchar();
																																for (i=0;i<clist[k].phased;i++) 
																																{
																																								if (slist[i] > 0 && h1[slist[i]] == '1') h1[slist[i]] = '0';
																																								else if (slist[i] > 0 && h1[slist[i]] == '0') h1[slist[i]] = '1';
																																}
																																clist[k].bestMEC = clist[k].MEC; clist[k].MEC =0; 
																																for (i=0;i<clist[k].frags;i++) 
																																{
																																								update_fragscore(Flist,clist[k].flist[i],h1); 
																																								clist[k].MEC += Flist[clist[k].flist[i]].currscore;
																																}
																																if (clist[k].MEC > clist[k].bestMEC) 
																																{
																																								for (i=0;i<clist[k].phased;i++)
																																								{
																																																if (slist[i] > 0 && h1[slist[i]] == '1') h1[slist[i]] = '0';
																																																else if (slist[i] > 0 && h1[slist[i]] == '0') h1[slist[i]] = '1';
																																								} clist[k].MEC = clist[k].bestMEC; 

																																}
																								} 
																								if (iter > 0 && clist[k].MEC > 0) fprintf(stdout,"component %d offset %d length %d phased %d  calls %d MEC %d cutvalue %f prevMEC %d bestMEC %d \n",k,clist[k].offset,clist[k].length,clist[k].phased,clist[k].calls,clist[k].MEC,cutvalue,clist[k].bestMEC,clist[k].lastMEC);
																								//if (iter > 0 && clist[k].MEC > 0) fprintf(stderr,"component %d offset %d length %d phased %d  calls %d MEC %d cutvalue %f prevMEC %d bestMEC %d \n",k,clist[k].offset,clist[k].length,clist[k].phased,clist[k].calls,clist[k].MEC,cutvalue,clist[k].bestMEC,clist[k].lastMEC);
																								if (iter >= 2 && clist[k].MEC > clist[k].lastMEC && iter%2==0 && clist[k].phased < 50 && clist[k].MEC > 0 ) 
																								{
																																for (i=0;i<clist[k].length;i++) 
																																{ 
																																								if (snpfrag[clist[k].offset+i].component == snpfrag[clist[k].offset].component) fprintf(stdout,"%c",h1[clist[k].offset+i]); else fprintf(stdout,"-"); 
																																} fprintf(stdout,"\n");
																																for (i=0;i<clist[k].length;i++) 
																																{ 
																																								if (snpfrag[clist[k].offset+i].component == snpfrag[clist[k].offset].component) fprintf(stdout,"%c",aaron[clist[k].offset+i]); else fprintf(stdout,"-"); 
																																} fprintf(stdout,"\n");
																																for (i=0;i<clist[k].length;i++) 
																																{ 
																																								if (snpfrag[clist[k].offset+i].component == snpfrag[clist[k].offset].component && slist[i] > 0) fprintf(stdout,"c"); else fprintf(stdout,"-"); 
																																} fprintf(stdout,"\n");
																								}
																}
								}
}




