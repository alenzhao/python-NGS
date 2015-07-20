/*
	*qsub chr22_qv-1.pl;
	perl -pi -e 's/error005/error01/g' chr22_error.pl; qsub chr22_error.pl
	perl -pi -e 's/error01/error02/g' chr22_error.pl ; qsub chr22_error.pl
	perl -pi -e 's/error02/error03/g' chr22_error.pl ; qsub chr22_error.pl
	perl -pi -e 's/error03/error05/g' chr22_error.pl ; qsub chr22_error.pl
	perl -pi -e 's/error05/error08/g' chr22_error.pl ; qsub chr22_error.pl
	perl -pi -e 's/error08/error1/g' chr22_error.pl ; qsub chr22_error.pl

	*
	*
	*/


#include "common.h"
int GIBBS =1;
int TREE =0; 
int AARON =1;
int STDERR =1;
int RANDOM_START=0;
int USE_HAP = 1;
int ERROR = 1;
int QV = 1;
int MINCUT =0;
int burnin = 20000;
int thinrate = 100;
int MCMCruns=0;
int SINGLE =0;
int UPDATE =1;
int UPDATETREE =1;
int EXAMPLE = 0; int ne = 20; int MF = 3;

int partitions(struct SNPfrags* snpfrag, int snps,int* slist, int N,int min,struct fragment* Flist,int fragments,FILE* mcout,int* nodes,char* hap);
int min_cut(struct SNPfrags* snpfrag, int snps,int* slist, int N,char* hap);
int tree_sampling(char* frags,char* sol,int fraction);
int sample_hap_rest(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1);


int compute_flist(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node, int mode)
{
								// mode ==1 corresponds to fragments linking two children and mode ==0 corresponds to fragments linking it to rest of SNPs 
								// mode ==1 sfrags   mode ==0 frags 
								int i=0,j=0,r=0,t=0,frags=0,one=0,zero=0; 
								int offset = clist[k].tree[node].offset, length = clist[k].tree[node].length,first=0, last=length-1;

								if (mode ==1)  // fragments linking two children case
								{
																for (i=0;i<length;i++) { if (clist[k].tree[node].bitvec[i] =='0') {last = i; if (first == 0) first = i; } }
								}

								for (i=snpfrag[first+offset].ff;i<fragments;i++)
								{
																// output fragments that connect both shores
																if (Flist[i].list[0].offset-offset  > last) break;
																if (Flist[i].component != snpfrag[offset].component) continue; one =0; zero=0;
																if (clist[k].tree[node].leaf =='1') 
																{
																								if(mode ==0) zero=1;
																								for (r=0;r<Flist[i].blocks;r++)
																								{
																																for (j=0;j<Flist[i].list[r].len;j++)
																																{
																																								t = Flist[i].list[r].offset+j - offset; 
																																								if (t < 0 || t >= length) continue;
																																								if (clist[k].tree[node].bitvec[t] == '1' && mode==0) one=1;  
																																								if (clist[k].tree[node].bitvec[t] == '1' && mode==1 && one ==0) one=1;  
																																								else if (clist[k].tree[node].bitvec[t] == '1' && mode==1 && one ==1) zero=1;  
																																}
																								}	
																}
																else	
																{
																								for (r=0;r<Flist[i].blocks;r++)
																								{
																																for (j=0;j<Flist[i].list[r].len;j++)
																																{
																																								t = Flist[i].list[r].offset+j - offset; 
																																								if (t < 0 || t >= length) { if (mode ==0) one = 1; continue; }
																																								if (clist[k].tree[node].bitvec[t] == '1' && mode==1) one=1;  
																																								if (clist[k].tree[node].bitvec[t] == '0' && mode ==1) zero=1; 
																																								if (clist[k].tree[node].bitvec[t] =='-' && mode ==0) one =1; else if (mode ==0) zero =1;
																																}
																								}
																} 	
																if (one ==1 && zero ==1 && mode ==0) { if (clist[k].tree[node].frags >0) clist[k].tree[node].flist[frags] = i; frags++; } 
																if (one ==1 && zero ==1 && mode ==1) { if (clist[k].tree[node].sfrags >0) clist[k].tree[node].sflist[frags] = i; frags++; } 
								} 
								if (clist[k].tree[node].frags ==-1 && mode ==0) clist[k].tree[node].frags = frags; 
								if (clist[k].tree[node].sfrags ==-1 && mode ==1) clist[k].tree[node].sfrags = frags; 
}

/* VERY IMPORTANT ERROR: BIT =0 or 1 different for sampleall and samplelocal *******/ 
/*************************************** GENERATE A SET OF PUTATIVE BLOCK UPDATES ****************************
	* general idea: strongly phased sets that are weakly connected/phased to rest of the snps 
	islands may be a good place to start 
	all global subsets of snps with low value of mincut to rest of snps 
	*/

int NJtree(struct SNPfrags* snpfrag, int snps,struct BLOCK* clist, int block,FILE* mcout,int* nodes)  
{
								// output each snp as a separate block--- identical to single SNP gibbs sampling 
								// clist[block].haplotype is set to 0 for all positions except '-' for snps from other components 
								fprintf(stderr," inside NJ tree code \n");
								int i=0,j=0,t=0,t1=0,t2=0,N= clist[block].phased;
								int pflag=0,W=0,k=0;
								int* slist = (int*)malloc(sizeof(int)*clist[block].phased);
								j=0; for (i=0;i<clist[block].length;i++) { if (clist[block].haplotype[i] != '-') { slist[j] = clist[block].offset+i; j++; } } 
								for (i=0;i<clist[block].length;i++) { if (clist[block].haplotype[i] != '-') clist[block].haplotype[i] = '0';  }
								for (i=0;i<N;i++)
								{
																j = slist[i];
																fprintf(mcout,"%d 1 1 1 \n",j);
																(*nodes)++;
								}

								i=0;	while(i<clist[block].length) 
								{
																if (clist[block].haplotype[i] == '-' || snpfrag[clist[block].offset+i].island =='0') i++; 
																else
																{
																								j=i; while (snpfrag[clist[block].offset+j].island =='1' && j < clist[block].length) j++;
																								fprintf(mcout,"%d %d %d ",clist[block].offset+i,j-i+1,j-i+1); 
																								for (t=i;t<=j;t++) fprintf(mcout,"1 "); fprintf(mcout,"\n");
																								(*nodes)++;
																								i = j+1;
																}
								}
								free(slist); 
								return 1;


								// add all cuts that partition the columns into two linear subsets 
								for (i=1;i<N-2;i++)
								{
																// compute the number of fragments connecting snp slist[i] and slist[i+1], if that is large then do not cut here
																// 

																fprintf(mcout,"%d %d %d ",slist[0],slist[i]-slist[0]+1,i+1);
																// print bitvector 
																t1=0;
																for (j=slist[0];j<= slist[i];j++)
																{
																								if (slist[t1] == j) { t1++; fprintf(mcout,"1 "); } else fprintf(mcout,"- "); 
																} fprintf(mcout,"\n");
																(*nodes)++;
								}
								free(slist); 
								return 1;

								// add cuts corresponding to well connected islands 
								// strategy for generating subsets online 
								// start from singletons and add subsets that are well connected and do not change phase over many iterations 
								//

								//								for (i=0;i<snps;i++) snpfrag[i].island = '0';
								//      for (i=0;i<fragments;i++) {for (j=0;j<Flist[i].blocks;j++){  for (k=0;k<Flist[i].list[j].len-1;k++) { snpfrag[Flist[i].list[j].offset+k].island = '1'; } } }
								//																fprintf(mcout,"%d 1 1 1 frags %d ",j,snpfrag[j].frags);				for(t=0;t<snpfrag[j].frags;t++) fprintf (mcout,"%d ",snpfrag[j].flist[t]); fprintf(mcout,"\n"); 



}
int compute_updatelist(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node)
{
								// find all other subsets (nodes) that are affected by flipping of this subset 
								int i=0,j=0,r=0,t=0,frags=0,one=0,zero=0,t1=0,t2=0,nodes=0; 
								int offset = clist[k].tree[node].offset, length = clist[k].tree[node].length; 

								for (t=0;t<clist[k].nodes;t++)
								{
																one =0;
																for (i=0;i<clist[k].tree[node].frags;i++)
																{		
																								j = clist[k].tree[node].flist[i]; 
																								for (t1=0;t1<clist[k].tree[t].frags;t1++) 
																								{
																																if (j == clist[k].tree[t].flist[t1]) { t1 = clist[k].tree[t].frags; i = clist[k].tree[node].frags; one =1;} 
																								}
																}
																if (one ==1) 
																{ 
																								if (clist[k].tree[node].updatenodes > 0 ) clist[k].tree[node].updatelist[nodes] = t; nodes++; 
																								continue;
																} 
																// if the subset being flipped 'node' contains some part of 't' and something else, then update 't'
																one =0; zero=0;
																for (i=0;i<clist[k].tree[node].length;i++)
																{
																								if (clist[k].tree[node].offset+i < clist[k].tree[t].offset || clist[k].tree[node].offset+i >= clist[k].tree[t].offset+clist[k].tree[t].length-1) zero=1;
																								else if (clist[k].tree[node].bitvec[i] != '-' && clist[k].tree[t].bitvec[clist[k].tree[node].offset+i-clist[k].tree[t].offset] != '-') one =1; 
																}
																if (one ==1 && zero == 1) 
																{ 
																								if (clist[k].tree[node].updatenodes > 0 ) clist[k].tree[node].updatelist[nodes] = t; nodes++; 
																} 

								} 
								clist[k].tree[node].updatenodes = nodes; return 1;
}



/************************************* CODE FOR DOINF TREE BASED EXACT SAMPLING **********************************************/

int compute_flipprob(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1,int update)
{
								int i=0,j=0,t=0, t1=0,t2=0,flip=0,bit=0;				double p[4], p11_00=0,p10_01=0,a,b;
								if (clist[k].tree[node].frags ==0) return 0;
								for (i=0;i<clist[k].tree[node].frags;i++)
								{
																j = clist[k].tree[node].flist[i]; p[0]=p[1]=p[2]=p[3]=0;
																for (t1=0;t1<Flist[j].blocks;t1++)
																{
																								for (t2=0;t2<Flist[j].list[t1].len;t2++)
																								{
																																t=Flist[j].list[t1].offset + t2; 
																																if (t-clist[k].tree[node].offset < 0 || t-clist[k].tree[node].offset >= clist[k].tree[node].length ) bit =0;
																																else  if (clist[k].tree[node].bitvec[t-clist[k].tree[node].offset] =='-') bit =0;
																																else bit =1; 
																																b = log10(1-Flist[j].list[t1].pv[t2]); a = log10(Flist[j].list[t1].pv[t2]);
																																if (PMODEL ==0) { b = log10(1-snpfrag[t].pv); a = log10(snpfrag[t].pv);}
																																if (bit ==1)
																																{
																																								if (h1[t] == Flist[j].list[t1].hap[t2]) { p[0] +=b; p[1] +=b; p[2] += a; p[3] +=a; }
																																								else { p[0] +=a; p[1] +=a; p[2] += b; p[3] +=b; }
																																}
																																else 
																																{
																																								if (h1[t] == Flist[j].list[t1].hap[t2]) { p[0] +=b; p[2] +=b; p[1] += a; p[3] +=a; }
																																								else { p[0] +=a; p[1] +=b; p[2] += a; p[3] +=b; }
																																}
																								}
																} 
																if (p[0] > p[3]) p11_00 += p[0] + log10(1+pow(10,p[3]-p[0])); else  p11_00 += p[3] + log10(1+pow(10,p[0]-p[3]));
																if (p[1] > p[2]) p10_01 += p[1] + log10(1+pow(10,p[2]-p[1])); else  p10_01 += p[2] + log10(1+pow(10,p[1]-p[2]));
								}
								if (update ==1) { clist[k].tree[node].p00 = p11_00; clist[k].tree[node].p01 = p10_01; } 
								if (update ==0) { clist[k].tree[node].P00 = p11_00; clist[k].tree[node].P01 = p10_01; } 

}



int sample_hap_rest(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1)
{
								// choose between flipping a subset of SNPs w.r.t remaining haplotype of clist[k]  
								int i=0,j=0,flip=0;	double a;		if (clist[k].tree[node].frags ==0) return 0;

								if (UPDATE ==0) compute_flipprob(snpfrag,snps,Flist,fragments,clist,k,node,h1,1); 
								else if (clist[k].tree[node].p00 ==100 && clist[k].tree[node].p01 ==100) compute_flipprob(snpfrag,snps,Flist,fragments,clist,k,node,h1,1);

								if (UPDATE == 3) // for checking if updating is working properly  
								{
																compute_flipprob(snpfrag,snps,Flist,fragments,clist,k,node,h1,0); 
																if (clist[k].tree[node].P00 > clist[k].tree[node].p00 + 0.02 || clist[k].tree[node].P00 < clist[k].tree[node].p00 - 0.02) { fprintf(stdout,"problem with probabilities %f %f \n",clist[k].tree[node].P00,clist[k].tree[node].p00); exit(0);}
																fprintf(stdout,"problem with probabilities %f %f \n",clist[k].tree[node].P00,clist[k].tree[node].p00); 
								}

								a = log10(drand48()); 	if (clist[k].tree[node].p01 >= a + clist[k].tree[node].p00) flip =1; else flip =0;		
								clist[k].tree[node].potflips++;
								if (flip ==1)
								{
																if (EXAMPLE && ne > 200) 
																{ 
																								fprintf(stdout,"flipping block %d offset %d end %d pr0 %f pr1 %f %f \n",node,clist[k].tree[node].offset,clist[k].tree[node].offset+clist[k].tree[node].length,clist[k].tree[node].p01,clist[k].tree[node].p00,a);
																								for (i=0;i<clist[k].length;i++) fprintf(stdout,"%c",h1[clist[k].offset+i]); fprintf(stdout,"\n");
																}
																//                compute_flipprob(snpfrag,snps,Flist,fragments,clist,k,node,h1); fprintf(stdout,"flipping block %d offset %d end %d pr0 %f pr1 %f %f \n",node,clist[k].tree[node].offset,clist[k].tree[node].offset+clist[k].tree[node].length,clist[k].tree[node].p01,clist[k].tree[node].p00,a); getchar();

																clist[k].tree[node].flips++;
																for (i=0;i<clist[k].tree[node].length;i++) 
																{
																								if (clist[k].tree[node].bitvec[i] !='-') h1[clist[k].tree[node].offset+i] = (char)(97 - (int)h1[clist[k].tree[node].offset+i]);
																} 
																if (UPDATE ==1) 
																{
																								for (i=0;i<clist[k].tree[node].updatenodes;i++) compute_flipprob(snpfrag,snps,Flist,fragments,clist,k,clist[k].tree[node].updatelist[i],h1,1);
																}
								}
								return flip;
}

void generate_example_2(int n)
{
								int i=0;
								char mfile[100]; sprintf(mfile,"example-%d.matrix",n);
								FILE* fp = fopen(mfile,"w");
								fprintf(fp,"%d %d\n",2*n+2,n);
								for(i=0;i<n/2;i++) fprintf(fp,"11543222 %d 00\n",i+1);
								for(i=0;i<n/2;i++) fprintf(fp,"11543222 %d 00\n",i+1);
								fprintf(fp,"12344444 %d 01\n",n/2);
								fprintf(fp,"12344444 %d 01\n",n/2);
								for(i=0;i<n/2-1;i++) fprintf(fp,"11543222 %d 00\n",n/2+i+1);
								for(i=0;i<n/2-1;i++) fprintf(fp,"11543222 %d 00\n",n/2+i+1); 
								fclose(fp);

								sprintf(mfile,"example-%d.phase",n);     fp = fopen(mfile,"w");
								fprintf(fp,"BLOCK: offset: 1 len: %d phased: %d\n",n,n);
								for(i=0;i<n/2;i++) fprintf(fp,"11543222 0\t 1\n");
								for(i=0;i<n/2;i++) fprintf(fp,"11543222 1\t 0\n");
								fclose(fp);

}

void generate_example(int n,int M)
{
								int i=0,j=0;
								char mfile[100]; sprintf(mfile,"example-%d.matrix",n);
								FILE* fp = fopen(mfile,"w");
								fprintf(fp,"%d %d\n",M*n,n);
								for (j=0;j<M;j++){         for(i=0;i<n/2;i++) fprintf(fp,"11543222 %d 00\n",i+1); } 
								for (j=0;j<M;j++) fprintf(fp,"12344444 %d 01\n",n/2);
								for (j=0;j<M;j++){        for(i=0;i<n/2-1;i++) fprintf(fp,"11543222 %d 00\n",n/2+i+1); } 
								fclose(fp);

								sprintf(mfile,"example-%d.phase",n);     fp = fopen(mfile,"w");
								fprintf(fp,"BLOCK: offset: 1 len: %d phased: %d\n",n,n);
								for(i=0;i<n/2;i++) fprintf(fp,"11543222 0\t 1\n");
								for(i=0;i<n/2;i++) fprintf(fp,"11543222 1\t 0\n");
								fclose(fp);

}

int main(int argc, char** argv)
{

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
								int* H0011 = (int*)malloc(4*snps);
								int* H0110 = (int*)malloc(4*snps);
								char* tree_besthap = (char*)malloc(snps+1);  char* besthap_mec = (char*)malloc(snps+1);
								for (i=0;i<snps;i++) { h1[i] = '-'; tree_hap[i]='-'; bn[i] = -1; aaron[i] = '-'; cn[i] = -1; H0011[i] = 0; H0110[i] =0;} // bn is component no, hold is old haplotype 
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
								fprintf(stderr,"input haplotype MEC %d calls %d log likelihood %f press enter to proceed\n",miscalls,calls,ll); //getchar();
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
																								//if (clist[k].phased > 300) continue;
								//																if (clist[k].phased > 50) continue;
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
																								if ((iter%100 ==0 && GIBBS ==1 && EXAMPLE ==0) || (iter%(snps*snps*((int)log2(snps))) ==0 && EXAMPLE ==1 && iter > 100000))
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
																								if (iter%2000==0 && iter > 0 && ERROR ==0 && SINGLE ==0) { fprintf(stdout,"iter %d OUTPUTTING BEST SOLUTION BY MEC USING ALL FRAGMENTS \n",iter); print_blocks(clist,components,tree_besthap,besthap_mec,h1,cn,Flist,fragments,bestscore_mec,snpfrag); exit(0); }
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
																								//								compare_Flist_hap(snpfrag,snps,Flist,fragments,aaron,Z,QV);
																								}
																								//Z++; update_posterior(Flist,fragments,h1,Z);
																								fprintf(stderr,"updating switch error estimates %d \n",Z);
																								switcherror_estimates(snpfrag,H0011,H0110, h1,clist,components,Z); Z++; 
																								
																}
								}
}


int tree_sampling_blockwise(char* frags,char* sol,int fraction)
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
								float errprob = (double)miscalls/(double)calls;  if (EXAMPLE) errprob = 0.05;

								if (STDERR) fprintf(stderr,"input haplotype MEC %d calls %d log likelihood %f\n",miscalls,calls,ll); 
								//								if (EXAMPLE) getchar();

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
								//	for (k=0;k<components;k++) { if (GIBBS ==1)   compute_tree_prob(snpfrag,snps,Flist,fragments,clist,k,0,h1,HAP,prob,ptemp0); } 

								float rd1 = drand48();
								int maxiter=200000;
								char mincut_temp[200]; sprintf(mincut_temp,"%s.mincut.output-%f",frags,rd1);
								int round = 0; int firstround =1;
								for (round =0; round < 2;round++) { 
																for (k=0;k<components;k++) // COMPUTATION OF TREE FOR EACH COMPONENT 
																{
																								if (clist[k].phased > 200 && EXAMPLE ==0 && round ==0) continue;
																								if (clist[k].phased <= 200 && EXAMPLE ==0 && round ==1) continue;
																								if (clist[k].phased < 50) maxiter = 50000; else maxiter = 200000;
																								if (clist[k].phased < 10) maxiter = 5000; 
																								if (clist[k].phased < 50 && SINGLE ==1) maxiter = 1000000; 
																								if (k != 12) continue; maxiter = 1000000;
																								fprintf(stdout,"MCMC for component %d length %d phased %d %d...%d initial MEC %d LL %f\n",k,clist[k].length,clist[k].phased,clist[k].offset,clist[k].offset+clist[k].length-1,clist[k].MEC,clist[k].LL); 
																								fprintf(stderr,"MCMC for component %d length %d phased %d %d...%d initial MEC %d LL %f\n",k,clist[k].length,clist[k].phased,clist[k].offset,clist[k].offset+clist[k].length-1,clist[k].MEC,clist[k].LL);
																								for (iter=0;iter<maxiter;iter++)
																								{
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
																																								if (clist[k].tree[t].phased >1 && RANDOM_START ==1 && firstround ==1) continue;

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
																																if ( (iter%1000==0 && iter > 0 && iter < 10000 ))// || iter > 10000 && iter%10000==0)
																																{ 
																																								if (SINGLE ==0 && UPDATETREE ==1) { clist[k].treecompute ='1'; clist[k].dealloc ='1'; clist[k].lastMEC = clist[k].bestMEC;	}
																																								if (SINGLE ==-1) 
																																								{
																																																for (i=0;i<clist[k].nodes;i++) { if (clist[k].tree[i].flips > 0) fprintf(stdout,"block %d offset %d phased %d potflips %d flips %d \n",i,clist[k].tree[i].offset,clist[k].tree[i].phased,clist[k].tree[i].potflips,clist[k].tree[i].flips); }
																																								} firstround =0; 
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
																																if ((iter%100 ==0 && GIBBS ==1 && EXAMPLE ==0) || (iter%100 ==0 && EXAMPLE ==1))
																																{
																																								if (STDERR) fprintf(stderr,"component %d offset %d phased %d length %d MCMC %6d %d %d %f BESTMEC %d Log-lh %f bestLL %f\n",k,clist[k].offset,clist[k].phased,clist[k].length,iter,clist[k].MEC,clist[k].calls,(double)clist[k].MEC/(double)clist[k].calls,clist[k].bestMEC,clist[k].LL,clist[k].bestLL);
																																								fprintf(stdout,"component %d offset %d phased %d length %d MCMC %6d %d %d %f BESTMEC %d Log-lh %f bestLL %f\n",k,clist[k].offset,clist[k].phased,clist[k].length,iter,clist[k].MEC,clist[k].calls,(double)clist[k].MEC/(double)clist[k].calls,clist[k].bestMEC,clist[k].LL,clist[k].bestLL);
																																								if (clist[k].bestMEC == MF && EXAMPLE ==1) { fprintf(stdout,"sampled best haplotype iter %d \n",iter); exit(0);}
																																								if (clist[k].bestMEC == 0 && iter >= 1000) maxiter = 1000;

																																}
																								}
																								fprintf(stdout,"iter %d OUTPUTTING BEST SOLUTION BY MEC USING ALL FRAGMENTS \n",iter); print_block(clist,k,aaron,besthap_mec,h1,cn,Flist,fragments,bestscore_mec,snpfrag);// exit(0);
																								fprintf(stdout,"MCMC for component %d length %d phased %d %d...%d final MEC %d LL %f\n",k,clist[k].length,clist[k].phased,clist[k].offset,clist[k].offset+clist[k].length-1,clist[k].bestMEC,clist[k].bestLL); 
																}
								}
}


/************************************* CODE FOR DOINF TREE BASED EXACT SAMPLING **********************************************
*/
int partitions(struct SNPfrags* snpfrag, int snps,int* slist, int N,int min,struct fragment* Flist,int fragments,FILE* mcout,int* nodes,char* hap)
{
								int i=0,j=0,k=0,t=0,one,zero,r; int* slist1; int* slist2; char* bitvec; int length = slist[N-1]-slist[0]+1, offset= slist[0], first = -1,last=0;
								bitvec = (char*)malloc(length+1); for (i=0;i<length;i++) bitvec[i] ='-';
								int frags=0;				int flag =1;
								//						fprintf(stdout,"offset %d length %d SIZE %d ",slist[0],length,N);
								if (flag) { 								fprintf(mcout,"%d %d %d ",slist[0],length,N); } 
								//								for (i=0;i<N;i++) fprintf(stdout,"%d ",slist[i]-slist[0]); fprintf(stdout,"\n");
								for (i=0;i<N;i++) bitvec[slist[i]-slist[0]] = '1';
								//for(i=0;i<length;i++) fprintf(stdout,"%c ",bitvec[i]); // fprintf(stdout," ");
								if (N > min) 
								{
																min_cut(snpfrag,snps,slist,N,hap);  // compute min-cut of graph represented by vertices in slist
																for (i=0;i<N;i++) { if (slist[i] <0) j++; } 
																if ( j < N/2) 
																{
																								for (i=0;i<N;i++) { if (slist[i] <0) bitvec[-1*slist[i]-1-offset] = '0'; } 
																}
																else
																{
																								for (i=0;i<N;i++) { if (slist[i] >=0 ) bitvec[slist[i]-offset] = '0'; }
																}
																if (flag) { for(i=0;i<length;i++) fprintf(mcout,"%c ",bitvec[i]);} // fprintf(stdout," first %d last %d\n",first,last);
																if (flag) fprintf(mcout,"\n");
																if (flag) *nodes = *nodes +1;
																free(bitvec);

																slist1 = (int*)malloc(4*j); k=0;
																for (i=0;i<N;i++) { if (slist[i] <0) {slist1[k] = -1*slist[i]-1; k++; }}
																partitions(snpfrag,snps,slist1,k,min,Flist, fragments,mcout,nodes,hap); 
																free(slist1); 

																slist2 = (int*)malloc(4*(N-j)); k=0;
																for (i=0;i<N;i++) { if (slist[i] >= 0) {slist2[k] = slist[i]; k++; }}
																partitions(snpfrag,snps,slist2,k,min,Flist,fragments,mcout,nodes,hap); 
																free(slist2); 
								} //getchar();
								else 
								{ 
																if (flag) { for(i=0;i<length;i++) fprintf(mcout,"%c ",bitvec[i]); }
																if (flag) fprintf(mcout,"\n");
																if (flag) *nodes = *nodes +1;
																free(bitvec);
								}
}
int phase_score(char* hap,int i, int j, char* p )
{
								if (USE_HAP ==0 && p[0] == p[1]) return 1; 
								else if (USE_HAP ==0 && p[0] != p[1]) return -1; 
								else if (USE_HAP ==0) return 1;

								if (hap[i] == hap[j] && p[0] == p[1] && USE_HAP ==1) return 1; 
								else if (hap[i] != hap[j] && p[0] != p[1] && USE_HAP==1) return 1; 
								else if (USE_HAP ==1) return -1; 
}

int min_cut(struct SNPfrags* snpfrag, int snps,int* slist, int N,char* hap)  
{
								// compute min-cut of graph represented by vertices in slist of size N 
								// this function returns a list of SUBSETS in order that represent consecutive MIN-CUTS
								// struct edge* telist; int tedges; // temporary edge list and number of edges for MIN CUT computation 
								int i=0,j=0,k=0,pflag=0; float W=0;
								double rd = drand48();
								if (pflag) { fprintf(stdout,"calling MIN CUT of size %d \n",N);	for (i=0;i<N;i++) fprintf(stdout,"%d %d ",i,slist[i]); fprintf(stdout,"\n");}
								// construct adjacency list representation of the smaller induced subgraph from BIG adjacency list 'elist'
								for (i=0;i<N;i++) snpfrag[slist[i]].node_id = rd;  // if node_id == rd it is part of graph
								for (i=0;i<N;i++)
								{
																snpfrag[slist[i]].tedges=0; k=-1;
																for (j=0;j<snpfrag[slist[i]].edges;j++) 
																{
																								if ( snpfrag[snpfrag[slist[i]].elist[j].snp].node_id != rd) continue;
																								if (k != snpfrag[slist[i]].elist[j].snp) { snpfrag[slist[i]].tedges++; k = snpfrag[slist[i]].elist[j].snp; } 
																}
																if (snpfrag[slist[i]].tedges ==1) 
																{
																								// simple min cut found here 
																								//																slist[i] = -1*slist[i]-1; return 1; 
																}
								}

								for (i=0;i<N;i++)
								{
																snpfrag[slist[i]].tedges=0; k=-1;
																for (j=0;j<snpfrag[slist[i]].edges;j++) 
																{
																								if ( snpfrag[snpfrag[slist[i]].elist[j].snp].node_id != rd) continue;
																								if (k != snpfrag[slist[i]].elist[j].snp) 
																								{ 
																																snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].snp = snpfrag[slist[i]].elist[j].snp; 
																																k = snpfrag[slist[i]].elist[j].snp; 
																																W = phase_score(hap,slist[i],k,snpfrag[slist[i]].elist[j].p);
																																snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].w = W; 
																																//                            snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].w = 1; 
																																snpfrag[slist[i]].tedges++;  
																								} 
																								else if (k == snpfrag[slist[i]].elist[j].snp) 
																								{
																																snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges-1].w += phase_score(hap,slist[i],k,snpfrag[slist[i]].elist[j].p); 
																																//																														snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges-1].w += 1;
																								}
																}
																if (pflag) { fprintf(stdout,"NODE %d ",slist[i]); for (j=0;j<snpfrag[slist[i]].tedges;j++) fprintf(stdout,"(%d,%d) ",snpfrag[slist[i]].telist[j].snp,snpfrag[slist[i]].telist[j].w); fprintf(stdout,"\n");}
								} // adjacency matrix is now captured in 'telist' and # of edges in 'tedges'
								// if hap is being used then make edges with negative weight to be 0, if not make all edge weights positive 
								for (i=0;i<N;i++)
								{
																//								    for (j=0;j<snpfrag[slist[i]].tedges;j++) { if (snpfrag[slist[i]].telist[j].w < 0 && USE_HAP ==1) snpfrag[slist[i]].telist[j].w = 0; } 
																for (j=0;j<snpfrag[slist[i]].tedges;j++) { if (snpfrag[slist[i]].telist[j].w < 0 && USE_HAP ==0) snpfrag[slist[i]].telist[j].w = -1*snpfrag[slist[i]].telist[j].w; } 
								}

								// NOW apply SIMPLE MIN CUT algorithm of wagner ESA 1994 to obtain 2 shores
								for (i=0;i<N;i++) snpfrag[slist[i]].parent = slist[i];  int V = N, Asize; int bscore,snp_add;
								//int* mergeorder = (int*)malloc(4*(2*N)); 
								int last1,last2, curr_cut,best_cut=100000,best_point;
								char* mincut = (char*)malloc(N); char* bestmincut = (char*)malloc(N); for (i=0;i<N;i++) bestmincut[i] = '0';
								int size_small,best_small=0,secondlast=0,last=0;
								while (V >1)
								{
																Asize = 1;
																for (i=0;i<N;i++) snpfrag[slist[i]].Aset = '0';  snpfrag[slist[0]].Aset = '1'; 
																//if (pflag) { for (i=0;i<N;i++) fprintf(stdout,"P(%d,%d) ",slist[i],snpfrag[slist[i]].parent); fprintf(stdout,"\n"); } 
																for (i=0;i<N;i++) mincut[i] = '0'; mincut[0] = '1';
																for (i=0;i<N;i++) snpfrag[slist[i]].score  = 0; bscore =-1000;
																for (j=0;j<snpfrag[slist[0]].tedges;j++)
																{
																								if (snpfrag[snpfrag[slist[0]].telist[j].snp].Aset =='0') snpfrag[snpfrag[snpfrag[slist[0]].telist[j].snp].parent].score += snpfrag[slist[0]].telist[j].w;
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
																								if (pflag) fprintf(stdout,"(%d %d)",slist[snp_add],bscore);
																								if (Asize == V-1) secondlast = snp_add;  if (Asize == V) last = snp_add;
																}
																/*
																			while (Asize < V)
																			{
																// find most tightly linked vertex to set Aset and add it to Aset 
																for (i=0;i<N;i++) snpfrag[slist[i]].score  = 0; bscore =-1000;
																for (i=0;i<N;i++) 
																{
																if (snpfrag[slist[i]].Aset  == '0') continue;
																for (j=0;j<snpfrag[slist[i]].tedges;j++)													
																{
																if (snpfrag[snpfrag[slist[i]].telist[j].snp].Aset =='0') snpfrag[snpfrag[snpfrag[slist[i]].telist[j].snp].parent].score += snpfrag[slist[i]].telist[j].w;
																}
																}
																for (i=0;i<N;i++) { if (snpfrag[slist[i]].score > bscore && snpfrag[slist[i]].Aset  =='0') {bscore = snpfrag[slist[i]].score; snp_add = i; } } 
																for (i=0;i<N;i++) if (snpfrag[slist[i]].parent == slist[snp_add]) snpfrag[slist[i]].Aset = '1';
																if (Asize < V-1) { for (i=0;i<N;i++) { if (snpfrag[slist[i]].parent == slist[snp_add] ) mincut[i] = '1'; }}
																Asize++; 
																if (pflag) fprintf(stdout,"(%d %d)",slist[snp_add],bscore);
																if (Asize == V-1) secondlast = snp_add;  if (Asize == V) last = snp_add;
																}
																*/
																// also maintain last two vertices added and they should be merged for next round of mincut computation 
																if (slist[secondlast] > slist[last]) snpfrag[slist[secondlast]].parent = snpfrag[slist[last]].parent; 
																else  snpfrag[slist[last]].parent = snpfrag[slist[secondlast]].parent;
																for (i=0;i<N;i++) snpfrag[slist[i]].parent =  snpfrag[snpfrag[slist[i]].parent].parent;
																curr_cut = bscore;  
																size_small =0; for (i=0;i<N;i++) {if (mincut[i] == '0') size_small++; } if (size_small >= N/2) size_small = N-size_small;
																for (i=0;i<N;i++) snpfrag[slist[i]].parent =  snpfrag[snpfrag[slist[i]].parent].parent;
																//				if (curr_cut < best_cut || (size_small > best_small && curr_cut*best_small < best_cut*size_small)) 
																//																if (curr_cut < best_cut || (size_small > best_small && curr_cut <= best_cut && GIBBS ==1)) 
																if (curr_cut < best_cut || (size_small > best_small && curr_cut <= best_cut )) 
																{ 
																								for (i=0;i<N;i++) bestmincut[i]= mincut[i];  best_small = size_small;   
																								best_cut = curr_cut; best_point = V; 
																}
																V--;
								} 
								//								fprintf(stderr,"slist %d MIN-CUT value %d point %d small %d N %d \n",slist[0],best_cut,best_point,best_small,N); 
								//								 fprintf(stderr,"slist %d MIN-CUT value %d point %d small %d N %d \n",slist[0],best_cut,best_point,best_small,N); 
								if (best_small ==0) 
								{ fprintf(stderr,"slist %d MIN-CUT value %d point %d small %d N %d \n",slist[0],best_cut,best_point,best_small,N); 

																{ fprintf(stderr,"BLOCK1 "); for (i=0;i<N;i++) { if (bestmincut[i] == '0') fprintf(stderr,"%d ",slist[i]);} fprintf(stderr,"\n");}
																{ fprintf(stderr,"BLOCK2 "); for (i=0;i<N;i++) { if (bestmincut[i] == '1') fprintf(stderr,"%d ",slist[i]);} fprintf(stderr,"\n");}
																getchar();
								}
								// compute the actual min cut 
								if (pflag) fprintf(stdout,"MIN-CUT value %d point %d \n",best_cut,best_point); 
								if (pflag) for (i=0;i<N;i++) fprintf(stdout,"%c",bestmincut[i]) ;
								if (pflag) { fprintf(stdout,"BLOCK1 "); for (i=0;i<N;i++) { if (bestmincut[i] == '0') fprintf(stdout,"%d ",slist[i]);} fprintf(stdout,"\n");}
								if (pflag) { fprintf(stdout,"BLOCK2 "); for (i=0;i<N;i++) { if (bestmincut[i] == '1') fprintf(stdout,"%d ",slist[i]);} fprintf(stdout,"\n");}
								for(i=0;i<N;i++) { if (bestmincut[i] == '1') slist[i]= -1*slist[i]-1; }

								free(mincut); free(bestmincut); return best_cut;
}
