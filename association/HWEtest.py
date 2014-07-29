#!/usr/bin/python2.4
import sys, os, glob, string, subprocess,time, math, re, compiler, random
from combinatorial import ncrbase2

# calculates HWE pvalue and returns #of expected hets, p-value for excess heterozygotes, and pvalue for excess homozygotes
# adapted from Abecasis paper in AJHG (exact tests of HWE)  
def HWE_exact(refs,hets,alts):
        N = refs + hets + alts; n0 = alts*2 + hets; n01 = hets;
        # N is sample size, n0 is number of minor allele, n01 is number of heterozygotes 
        # Plow = P(N01 <= n01 | N,n0)
        if n0 > N: n0 = 2*N-n0;
        p = float(n0)/(2*N);    plow = 0; phigh =0; constant = ncrbase2(2*N,n0);
        if n0%2 ==1: # odd
                for N01  in range(1,n01+1,2):
                        N0 = n0; N00 = (n0-N01)/2; N11 = N - N01 - N00; N1 = 2*N-N0;
                        plow += math.pow(2,N01 + ncrbase2(N,N00) + ncrbase2(N-N00,N01) - constant);
                for N01  in range(n01,n0+1,2):
                        N0 = n0; N00 = (n0-N01)/2; N11 = N - N01 - N00; N1 = 2*N-N0;
                        phigh += math.pow(2,N01 + ncrbase2(N,N00) + ncrbase2(N-N00,N01) - constant);

        if n0%2 ==0: # even
                for N01  in range(0,n01+1,2):
                        N0 = n0; N00 = (n0-N01)/2; N11 = N - N01 - N00; N1 = 2*N-N0;
                        plow += math.pow(2,N01 + ncrbase2(N,N00) + ncrbase2(N-N00,N01) - constant);
                for N01  in range(n01,n0+1,2):
                        N0 = n0; N00 = (n0-N01)/2; N11 = N - N01 - N00; N1 = 2*N-N0;
                        phigh += math.pow(2,N01 + ncrbase2(N,N00) + ncrbase2(N-N00,N01) - constant);
#       print 2*p*(1-p)*N,plow,phigh;
        return [2*p*(1-p)*N,plow,phigh];

print HWE_exact(100,0,80);
~                                                                                                                                                                                              
~                                
