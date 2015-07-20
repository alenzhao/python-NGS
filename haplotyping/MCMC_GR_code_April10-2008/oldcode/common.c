#include "common.h"

int compare_haps(struct BLOCK* clist, int components, char* orig, char* h1,struct SNPfrags* snpfrag, int snps)
{
								// if the two haplotypes differ at SNPs covered by single fragment, then we should not count these differences since there is no way we can reverse these errors
								int i=0,j=0,k=0,match=0,mismatch=0,flips=0,error=0,diffs=0,totaldiffs=0,matches=0;
								for (k=0;k<components;k++)
								{
																match=-1,mismatch=-1; diffs=0; matches=0;
																for (i=0;i<clist[k].length;i++)
																{
																								if (clist[k].haplotype[i] =='-') continue;
//																								if (snpfrag[clist[k].offset+i].frags <= 1) continue;
																								if (orig[clist[k].offset+i] == h1[clist[k].offset+i] && match ==-1 ) { matches++; match =1; }
																								else if (orig[clist[k].offset+i] == h1[clist[k].offset+i] && mismatch == 1 ) { mismatch =-1; matches++; match =1; } 
																								else if (orig[clist[k].offset+i] != h1[clist[k].offset+i] && match == 1 ) { match =-1; diffs++; mismatch = 1; } 
																								else if (orig[clist[k].offset+i] != h1[clist[k].offset+i] && mismatch ==-1 ) { diffs++; mismatch =1; }
																}
																if (diffs < matches) totaldiffs += diffs; else totaldiffs += matches;
								}
								return totaldiffs;
}

int compare_Flist_hap(struct SNPfrags* snpfrag, int snps, struct fragment* Flist,int fragments,char* h,int Z,int QV)
{
								int f=0,j=0,k=0,good=0,bad=0,flips=0,half=0,all=0,all_half=0,neg=0,goodcalls=0,badcalls=0,singletons=0;
								float good_ll=0,bad_ll=0,q=0;
								for (f=0;f<fragments;f++)
								{
																good=bad=0;
																if (Flist[f].blocks ==1 && Flist[f].list[0].len ==1) { singletons++; continue;}
																for (j=0;j<Flist[f].blocks;j++)
																{
																								for (k=0;k<Flist[f].list[j].len;k++)
																								{
																																if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good++; else bad++;
																								}
																}
																fprintf(stdout,"%d %s ",Flist[f].blocks,Flist[f].id);
																for (j=0;j<Flist[f].blocks;j++)
																{
																								fprintf(stdout,"%d %s ",Flist[f].list[j].offset,Flist[f].list[j].hap);
																} fprintf(stdout,"\n");
																for (j=0;j<Flist[f].blocks;j++)
																{
																								if (Z ==0)
																								{
																																for (k=0;k<Flist[f].list[j].len;k++)
																																{
																																								if ( (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k] && good < bad) || (h[Flist[f].list[j].offset+k] != Flist[f].list[j].hap[k] && good >= bad)) { fprintf(stdout,"%f ",Flist[f].list[j].pv[k]); badcalls++; bad_ll += log10(Flist[f].list[j].pv[k]); } 
																																								else { goodcalls++; good_ll += log10(Flist[f].list[j].pv[k]); } 
																																}
																								}
																								else
																								{
																																for (k=0;k<Flist[f].list[j].len;k++)
																																{
																																				if (QV != -1) { q = (int)(Flist[f].list[j].qv[k])-48; q /=10; q *= -1; q = pow(10,q); } 
                                    else q = Flist[f].list[j].pv[k]; 
																																				if (Flist[f].list[j].qv[k] =='1') fprintf(stdout,"BC %d %d %f %f %d FLIP %c\n",f,Flist[f].list[j].offset+k,q,Flist[f].list[j].post[k]/Z,snpfrag[Flist[f].list[j].offset+k].frags,Flist[f].list[j].qv[k]);
																																				else fprintf(stdout,"BC %d %d %f %f %d NORM %c\n",f,Flist[f].list[j].offset+k,q,Flist[f].list[j].post[k]/Z,snpfrag[Flist[f].list[j].offset+k].frags,Flist[f].list[j].qv[k]);
																																}
																								}

																}
								}
								if (Z ==0) {
								fprintf(stdout,"goodcalls %d good_LL %f avgLL %f     badcalls %d bad_LL %f avgLL %f \n",goodcalls,good_ll,good_ll/goodcalls,badcalls,bad_ll,bad_ll/badcalls);
								fprintf(stderr,"goodcalls %d good_LL %f avgLL %f     badcalls %d bad_LL %f avgLL %f \n",goodcalls,good_ll,good_ll/goodcalls,badcalls,bad_ll,bad_ll/badcalls);
								fprintf(stderr,"%d \n",singletons);
							 }
}

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

int mutate_Flist(struct fragment* Flist,int fragments,double errrate)
{
								int j=0,k=0,f=0,flips=0;
								for (f=0;f<fragments;f++)
								{
																for (j=0;j<Flist[f].blocks;j++)
																{
																								for (k=0;k<Flist[f].list[j].len;k++)
																								{
																																Flist[f].list[j].qv[k] = '0';
																																if (drand48() < errrate) 
																																{ 
																																								Flist[f].list[j].hap[k] = (char)(97- (int)Flist[f].list[j].hap[k]); flips++; 
																																								Flist[f].list[j].qv[k] = '1'; 
																																				//fprintf(stderr,"f %d POS %d %c\n",f,Flist[f].list[j].offset+k,Flist[f].list[j].qv[k]);
																																				//fprintf(stdout,"f %d POS %d %c\n",f,Flist[f].list[j].offset+k,Flist[f].list[j].qv[k]);
																																}
																								}
																}
								}
								return flips;
}

int correct_fragment(struct fragment* Flist,int f, char* h)
{
								int j=0,k=0,good=0,bad=0;
								for (j=0;j<Flist[f].blocks;j++)
								{
																for (k=0;k<Flist[f].list[j].len;k++)
																{
																								if (h[Flist[f].list[j].offset+k] == '-') continue;
																								if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good++; else bad++;
																}
								}
								if (good ==0 || bad ==0) return 1;
								for (j=0;j<Flist[f].blocks;j++)
								{
																for (k=0;k<Flist[f].list[j].len;k++)
																{
																								if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k] && good < bad) Flist[f].list[j].hap[k] = (char)(97- (int)Flist[f].list[j].hap[k]);
																								if (h[Flist[f].list[j].offset+k] != Flist[f].list[j].hap[k] && good >= bad) Flist[f].list[j].hap[k] = (char)(97- (int)Flist[f].list[j].hap[k]);
																}
								}
}


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
																								if (h[Flist[f].list[j].offset+k] == '-') continue;// { fprintf(stdout,"fragment error"); continue;}
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


int gibbs_sampling(char* frags,char* sol,int fraction);
int sample_block(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int start, int end,char* h1,int snps);


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
																																if (Flist[i].floppy ==1) continue;
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
								int j=0; for (j=0;j<blocks;j++) 
								{
												if (blist[j].length > 120 || blist[j].length < 40) continue; 
												print_block(blist,j,aaron,h1,current,bn,Flist,fragments,mecscore,snpfrag);
								}
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
																//if ( ( diffs == 0 || pos == diffs) ) return 1;

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

