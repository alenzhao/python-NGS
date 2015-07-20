

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#include "common.h"
#include "mcmcsampling.h"


int switcherror_estimates(struct SNPfrags* snpfrag, int* H0011, int* H0110,char* h1, struct BLOCK* clist, int components,int Z)
{
				int i=0,j=0,first,k=0;
				if (Z%500 ==0) { fprintf(stderr,"outputting switch errors \n"); } 
				for (k=0;k<components;k++)			
				{
								if (clist[k].phased > 50 ) continue;
								first = clist[k].offset;
								for (i=1;i<clist[k].length;i++)
								{
												if (snpfrag[clist[k].offset+i].component == snpfrag[clist[k].offset].component)
												{
																if (h1[first] == '0' && h1[clist[k].offset+i] =='0') H0011[first]++; 
																if (h1[first] == '1' && h1[clist[k].offset+i] =='1') H0011[first]++; 
																if (h1[first] == '1' && h1[clist[k].offset+i] =='0') H0110[first]++; 
																if (h1[first] == '0' && h1[clist[k].offset+i] =='1') H0110[first]++; 
																if (Z%500==0 && Z > 0) fprintf(stdout,"SER Z %d %d %d -- 00 %f 11 %f\n",Z,first,clist[k].offset+i,(float)H0011[first]/(float)(Z+1),(float)H0110[first]/(float)(Z+1));
																first = clist[k].offset+i;
												}
								}
								
				}
				if (Z%500==0 && Z > 0 ) { fprintf(stderr,"outputting switch errors \n"); }//getchar(); } 
}

int update_posterior(struct fragment* Flist,int fragments,char* h,int Z)
{
								int f=0,j=0,k=0,good=0,bad=0,flips=0,half=0,all=0,all_half=0,neg=0,zi=0;
								float p0=0,p1=0,sum=0;
								for (f=0;f<fragments;f++)
								{
																good=bad=0; p0=0,p1=0;zi=0;
																for (j=0;j<Flist[f].blocks;j++)
																{
																								for (k=0;k<Flist[f].list[j].len;k++)
																								{
																																if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good++; else bad++;
																																if (Z==1) Flist[f].list[j].post[k] =0;
																																if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) 
																																{ p0 += log10(1-Flist[f].list[j].pv[k]); p1 += log10(Flist[f].list[j].pv[k]); }
																																else
																																{ p1 += log10(1-Flist[f].list[j].pv[k]); p0 += log10(Flist[f].list[j].pv[k]); }
																																if (Z ==0) Flist[f].list[j].post[k] = 0;
																								}
																}
																if (p0 > p1) { p1 -= p0; p0= 0;} else { p0 -= p1; p1=0; } sum = pow(10,p0) + pow(10,p1);
																if (pow(10,p0) >= drand48()*sum) zi=0; else zi=1;
																for (j=0;j<Flist[f].blocks;j++)
																{
																								for (k=0;k<Flist[f].list[j].len;k++)
																								{
//																																if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) Flist[f].list[j].post[k] += p0-p1; else Flist[f].list[j].post[k] += p1-p0; 
																																//if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k] && zi ==0 ) Flist[f].list[j].post[k]++; 
																																//if (h[Flist[f].list[j].offset+k] != Flist[f].list[j].hap[k] && zi ==1 ) Flist[f].list[j].post[k]++; 
																																if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k] && zi ==1 ) Flist[f].list[j].post[k]--; 
																																if (h[Flist[f].list[j].offset+k] != Flist[f].list[j].hap[k] && zi ==0 ) Flist[f].list[j].post[k]--; 
																								}
																}	
								}
}

int sample_pair_hap(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1,int compute)
{
								// clist[k].bitvec tells which snp is in which block 000-1-000     111
								int i=0,j=0,t=0, t1=0,t2=0;				double p[4], p11_00=0,p10_01=0,a,b;

								for (i=0;i<clist[k].tree[node].sfrags;i++)
								{
																j = clist[k].tree[node].sflist[i]; p[0]=p[1]=p[2]=p[3]=0;
																//										fprintf(stdout,"\n%d %s ",Flist[j].blocks,Flist[j].id); 
																//for (t=0;t<Flist[j].blocks;t++) fprintf(stdout,"| %d %s %d | ",Flist[j].list[t].offset,Flist[j].list[t].hap,Flist[j].list[t].len); 
																for (t1=0;t1<Flist[j].blocks;t1++)
																{
																								for (t2=0;t2<Flist[j].list[t1].len;t2++)
																								{
																																t=Flist[j].list[t1].offset + t2; 
																																// clist[k].tree[i].bitvec = (char*)malloc(clist[k].tree[i].length+1);
																																if (t-clist[k].tree[node].offset < 0) continue;
																																if (t-clist[k].tree[node].offset >= clist[k].tree[node].length ) continue; //{ t1 = Flist[j].blocks; break; } 
																																b = log10(1-Flist[j].list[t1].pv[t2]); a = log10(Flist[j].list[t1].pv[t2]);
																																if (PMODEL ==0) { b = log10(1-snpfrag[t].pv); a = log10(snpfrag[t].pv);}
																																//																												fprintf(stdout,"%c %c %d |",Flist[j].list[t1].hap[t2],clist[k].tree[node].bitvec[t-clist[k].tree[node].offset],t); 
																																if (clist[k].tree[node].bitvec[t-clist[k].tree[node].offset] =='0') 
																																{
																																								if (h1[t] == Flist[j].list[t1].hap[t2]) { p[0] +=b; p[1] +=b; p[2] += a; p[3] +=a; }
																																								else { p[0] +=a; p[1] +=a; p[2] += b; p[3] +=b; }
																																}
																																else if (clist[k].tree[node].bitvec[t-clist[k].tree[node].offset] =='1') 
																																{
																																								if (h1[t] == Flist[j].list[t1].hap[t2]) { p[0] +=b; p[2] +=b; p[1] += a; p[3] +=a; }
																																								else { p[0] +=a; p[1] +=b; p[2] += a; p[3] +=b; }
																																}
																								}
																} 
																if (p[0] > p[3]) p11_00 += p[0] + log10(1+pow(10,p[3]-p[0])); else  p11_00 += p[3] + log10(1+pow(10,p[0]-p[3]));
																if (p[1] > p[2]) p10_01 += p[1] + log10(1+pow(10,p[2]-p[1])); else  p10_01 += p[2] + log10(1+pow(10,p[1]-p[2]));
																//														fprintf(stdout,"00 %f 01 %f 10 %f 11 %f || %f %f ",p[0],p[1],p[2],p[3],p11_00,p10_01);
								}
								//								fprintf(stdout,"\ncalling pairwise sampling 00 %f 10 %f\n %d %d \n",p11_00,p10_01,clist[k].tree[node].offset,clist[k].tree[node].length); 
								if (p11_00 > p10_01) { p10_01 -= p11_00; p11_00 = 0; } else { p11_00 -= p10_01; p10_01 = 0; } 
								if (compute ==1)
								{
																clist[k].tree[node].prob = clist[k].tree[clist[k].tree[node].ch1].prob + clist[k].tree[clist[k].tree[node].ch2].prob;
																if (p11_00 < p10_01) clist[k].tree[node].prob += p11_00 - p10_01; 
																//											fprintf(stderr," node %d children %f node %f %f total %f \n",node,clist[k].tree[clist[k].tree[node].ch1].prob + clist[k].tree[clist[k].tree[node].ch2].prob,p11_00,p10_01,clist[k].tree[node].prob); getchar();
																return 1;
																a = pow(10,p11_00) + pow(10,p10_01);
																clist[k].tree[node].prob += p11_00 - log10(a); 
								}
								b = pow(10,p10_01); a = drand48()*(pow(10,p11_00) + b); 
								if ( (a <= b && compute ==0) || (compute ==2 && p10_01 > p11_00)) 
								{
																// flip the '0' bits of bitvector 
																for (i=0;i<clist[k].tree[node].length;i++) 
																{
																								if (clist[k].tree[node].bitvec[i] =='0') h1[clist[k].tree[node].offset+i] = (char)(97 - (int)h1[clist[k].tree[node].offset+i]);
																}
								}
								return 1;
}

int sample_local_hap(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k, int node,char* h1,char** HAP,float* prob,float* probtemp,int compute)
								// compute is 1 we want to compute the relative prob of the haplotype 
								// compute is 0 we want to sample a local haplotype 
{
								int i=0,j=0,t=0, t1=0,t2=0,h=0,start = clist[k].tree[node].offset, last = clist[k].tree[node].length+start-1;
								//		if (phased < clist[k].tree[node].length) { fprintf(stderr," need to change code ... \n"); } 
								if (clist[k].tree[node].phased > 10) return -1;
								if (clist[k].tree[node].phased < 2 && compute ==0)		{ h1[clist[k].tree[node].offset] = '0'; return 0; } 
								if (clist[k].tree[node].phased < 2 && compute ==1)		{ clist[k].tree[node].prob = 0; return 0; }
								//fprintf(stderr,"iprob %f ",clist[k].tree[node].prob); 

								short* map = (short*)malloc(sizeof(short)*clist[k].tree[node].length);
								j=0; for (i=0;i<clist[k].tree[node].length;i++) { map[i] = -1; if (clist[k].tree[node].bitvec[i] != '-') { map[i] = j; j++; } } 
								//								j=0; for (i=0;i<clist[k].tree[node].length;i++) fprintf(stdout,"%d %d ",i,map[i]);
								int	haps = pow(2,clist[k].tree[node].phased); double a,b,max;
								for (h=0;h<haps;h++) prob[h] =0;

								//fprintf(stdout," local phasing .... %d\n",clist[k].tree[node].frags);
								for (i=0;i<clist[k].tree[node].sfrags;i++)
								{
																j = clist[k].tree[node].sflist[i]; for (h=0;h<haps;h++) probtemp[h]=0;
																//					fprintf(stdout,"\n%d %s ",Flist[j].blocks,Flist[j].id); for (t=0;t<Flist[j].blocks;t++) fprintf(stdout,"| %d %s %d | ",Flist[j].list[t].offset,Flist[j].list[t].hap,Flist[j].list[t].len); 
																for (t1=0;t1<Flist[j].blocks;t1++)
																{
																								for (t2=0;t2<Flist[j].list[t1].len;t2++)
																								{
																																t=Flist[j].list[t1].offset + t2; 
																																if (t-start < 0 || t > last) continue; 
																																if (map[t-start] == -1) continue;
																																b = log10(1-Flist[j].list[t1].pv[t2]); a = log10(Flist[j].list[t1].pv[t2]);
																																if (PMODEL ==0) { b = log10(1-snpfrag[t].pv); a = log10(snpfrag[t].pv); }
																																//																				fprintf(stdout,"M %d %d ",t-start,map[t-start]); fflush(stdout);
																																for (h=0;h<haps;h++) { if (HAP[h][map[t-start]] == Flist[j].list[t1].hap[t2]) probtemp[h] +=b; else probtemp[h] +=a; } 
																								}
																}
																for (h=0;h<haps/2;h++) 
																{ 
																								if (probtemp[h] > probtemp[haps-h-1])				probtemp[h] += log10(1+pow(10,probtemp[haps-h-1]-probtemp[h])); 
																								else			probtemp[h] = probtemp[haps-h-1] + log10(1+pow(10,probtemp[h]-probtemp[haps-h-1])); 
																}
																for (h=0;h<haps/2;h++) prob[h] += probtemp[h];
								}
								max = -100; for (i=0;i<haps/2;i++) { if (prob[i] > max) max = prob[i]; }
								b=0; for (i=0;i<haps/2;i++) b += pow(10,prob[i]-max);  

								if (compute ==2) // take the most likely haplotype 
								{
																for (i=0;i<haps/2;i++) { if (prob[i] >= max) h = i; }
																j=0; for (i=start;i<=last;i++) { if (clist[k].tree[node].bitvec[i-start] != '-') {h1[i] = HAP[h][j]; j++;} }
								}
								if (compute ==0)
								{
																a = drand48()*b;
																for (i=0;i<haps/2;i++) {  b = pow(10,prob[i]-max); if (a < b) { h = i; i= haps; } a -= b; }
																j=0; for (i=start;i<=last;i++) { if (clist[k].tree[node].bitvec[i-start] != '-') {h1[i] = HAP[h][j]; j++;} } 
																clist[k].tree[node].H = h; clist[k].tree[node].prob = prob[h]-max-log10(b); 
								}
								else if (compute ==1) // compute the prob of haplotype h1'''' 
								{
																j=0; for (i=last;i>= start;i--) { if (clist[k].tree[node].bitvec[i-start] != '-') j = 2*j + (int)(h1[i]-48); }
																if (j >= haps/2) j = haps-j-1; clist[k].tree[node].prob = prob[j]-max; // - log10(b); 
																//					fprintf(stderr,"cprob %f \n",clist[k].tree[node].prob);
								}
								free(map);

								//fprintf(stdout,"\n start %d end %d phased %d \n",start,last,phased);												   for (h=0;h<haps/2;h++) fprintf(stdout,"%s %f \n",HAP[h],prob[h]);														fprintf(stdout,"calling local sampling.... %f %f \n",max,max2); //getchar();

								//free(prob); free(probtemp); for (i=0;i<haps;i++) free(HAP[i]); free(HAP); 
								//for (i=0;i<haps;i++) fprintf(stdout,"%d %s\n",i,HAP[i]); 

}
int compute_tree_prob(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp)
{
								if (clist[k].tree[node].leaf =='1')		sample_local_hap(snpfrag,snps,Flist,fragments,clist,k,node,h1,HAP,prob,probtemp,1); 
								else 
								{
																compute_tree_prob(snpfrag,snps,Flist,fragments,clist,k,clist[k].tree[node].ch1,h1,HAP,prob,probtemp); 
																compute_tree_prob(snpfrag,snps,Flist,fragments,clist,k,clist[k].tree[node].ch2,h1,HAP,prob,probtemp); 
																sample_pair_hap(snpfrag,snps,Flist,fragments,clist,k,node,h1,1);
								}
}

int update_tree_prob(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp)
{
								int p = node;
								//								fprintf(stderr,"\n updating node probs..... \n");
								while (p >= 0)	
								{
																if (clist[k].tree[p].leaf =='1')		sample_local_hap(snpfrag,snps,Flist,fragments,clist,k,p,h1,HAP,prob,probtemp,1);
																else sample_pair_hap(snpfrag,snps,Flist,fragments,clist,k,p,h1,1); 
																//									fprintf(stderr,"node %d phased %d prob %f parent %d\n",p,clist[k].tree[p].phased,clist[k].tree[p].prob,clist[k].tree[p].parent);
																p = clist[k].tree[p].parent;
								} //getchar();
}

int sample_tree(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments,struct BLOCK* clist,int k,int node,char* h1,char** HAP,float* prob,float* probtemp,int ML) // ML ==2 iimplies maximum likelihood else sample 
{
								// sample haplotype for block k using the tree clist[k].tree 
								//fprintf(stdout,"node %d offset %d phased %d child1 %d child2 %d frags %d\n",node,clist[k].tree[node].offset,clist[k].tree[node].phased,clist[k].tree[node].ch1,clist[k].tree[node].ch2,clist[k].tree[node].frags);
								if (clist[k].tree[node].leaf =='1')		{ sample_local_hap(snpfrag,snps,Flist,fragments,clist,k,node,h1,HAP,prob,probtemp,ML); return 0;}
								sample_tree(snpfrag,snps,Flist,fragments,clist,k,clist[k].tree[node].ch1,h1,HAP,prob,probtemp,ML); 
								sample_tree(snpfrag,snps,Flist,fragments,clist,k,clist[k].tree[node].ch2,h1,HAP,prob,probtemp,ML); 
								sample_pair_hap(snpfrag,snps,Flist,fragments,clist,k,node,h1,ML);
}

int construct_tree(struct BLOCK* clist,int k,int parent, int st) // PREORDER TRAVERSAL ---> BINARY TREE 
{
								if (st >= clist[k].nodes) return -1; 
								else if ( clist[k].tree[st].leaf =='1') 
								{
																clist[k].tree[st].parent = parent; 
																if (clist[k].tree[parent].ch1 == -1) { clist[k].tree[parent].ch1 = st; construct_tree(clist,k,parent,st+1);}
																else if (clist[k].tree[parent].ch2 == -1) { clist[k].tree[parent].ch2 = st; construct_tree(clist,k,clist[k].tree[parent].parent,st+1); } 
																else construct_tree(clist,k,clist[k].tree[parent].parent,st);
								}
								else if ( clist[k].tree[parent].ch1 == -1) 
								{  
																clist[k].tree[parent].ch1 = st; clist[k].tree[st].parent = parent; 
																if (clist[k].tree[st].leaf =='1') construct_tree(clist,k,parent,st+1); else construct_tree(clist,k,st,st+1);
								} 
								else if ( clist[k].tree[parent].ch2 == -1) 
								{  
																clist[k].tree[parent].ch2 = st; clist[k].tree[st].parent = parent; 
																if (clist[k].tree[st].leaf =='1') construct_tree(clist,k,clist[k].tree[parent].parent,st+1); else construct_tree(clist,k,st,st+1);
								} 
								else construct_tree(clist,k,clist[k].tree[parent].parent,st);
}


/************************************* CODE FOR DOINF TREE BASED EXACT SAMPLING **********************************************/


int sample_haplotype(struct SNPfrags* snpfrag,int snps,struct fragment* Flist,int fragments, struct BLOCK* clist, int k,char* h1,char** HAP,float* prob,float* probtemp)
{
								int i=0,j=0,t=0, t1=0,t2=0,h=0,start = clist[k].offset, last = clist[k].length+start-1;
								int	haps = pow(2,clist[k].phased); 
								short* map = (short*)malloc(sizeof(short)*clist[k].length);
								j=0; for (i=0;i<clist[k].length;i++) { map[i] = -1; if (clist[k].haplotype[i] != '-') { map[i] = j; j++; } } 
								double a,b,max,max2;

								for (h=0;h<haps;h++) prob[h] =0;

								for (i=0;i<clist[k].frags;i++)
								{
																j = clist[k].flist[i]; for (h=0;h<haps;h++) probtemp[h]=0;
																for (t1=0;t1<Flist[j].blocks;t1++)
																{
																								for (t2=0;t2<Flist[j].list[t1].len;t2++)
																								{
																																t=Flist[j].list[t1].offset + t2; 
																																if (t-start < 0) continue; if ( t > last) continue; //{ t1 = Flist[j].blocks; break; } 
																																if (map[t-start] == -1) continue;
																																b = log10(1-Flist[j].list[t1].pv[t2]); a = log10(Flist[j].list[t1].pv[t2]);
																																if (PMODEL ==0) { b = log10(1-snpfrag[t].pv); a = log10(snpfrag[t].pv); }
																																//																				fprintf(stdout,"M %d %d ",t-start,map[t-start]); fflush(stdout);
																																for (h=0;h<haps;h++) { if (HAP[h][map[t-start]] == Flist[j].list[t1].hap[t2]) probtemp[h] +=b; else probtemp[h] +=a; } 
																								}
																}
																for (h=0;h<haps/2;h++) 
																{ 
																								if (probtemp[h] > probtemp[haps-h-1])				probtemp[h] += log10(1+pow(10,probtemp[haps-h-1]-probtemp[h])); 
																								else			probtemp[h] = probtemp[haps-h-1] + log10(1+pow(10,probtemp[h]-probtemp[haps-h-1])); 
																}
																for (h=0;h<haps/2;h++) prob[h] += probtemp[h];
								}
								max = -100; for (i=0;i<haps/2;i++) { if (prob[i] > max) max = prob[i];}
								b=0; for (i=0;i<haps/2;i++) b += pow(10,prob[i]-max); 
								//for (i=0;i<haps/2;i++) {	t = (int)(10000*pow(10,prob[i]-max)/b); if (t ==0) continue; for (j=0;j<clist[k].phased;j++) fprintf(stdout,"%c",HAP[i][j]); fprintf(stdout," | %d LL %f count %d \n",i,prob[i]-max,t);	} exit(0);

								a = drand48()*b;
								for (i=0;i<haps/2;i++) {  b = pow(10,prob[i]-max); if (a <= b) { h = i; i= haps; } a -= b; }
								j=0; for (i=start;i<=last;i++) { if (clist[k].haplotype[i-start] != '-') {h1[i] = HAP[h][j]; j++;} } 
								//free(prob); free(probtemp); for (i=0;i<haps;i++) free(HAP[i]); free(HAP); 
								free(map);

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

