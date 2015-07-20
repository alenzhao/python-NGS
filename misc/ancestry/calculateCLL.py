#! /usr/bin/env python

import sys, os, glob, string, subprocess,time, math, re, compiler, random

## input files are allelefreqmatrix for K populations (S x K), plus allele information for each SNP 0/1
## admixture coefficients for I individuals (I x K) 
## and genotypes in plink ped format for I individuals (I x S) 

## calculate CLL for each individual using given admixture proportions


##python calculate_CLL.py 11pops.estimate.AF proportions test.249075.ped > r

## read allele frequency matrix for K populations and allele per SNP 
def read_matrix(affile):
        afmatrix = []; alleles = [];
        File = open(affile,'r');
        for line in File:
                if line[0]=='#':continue;
                snp = line.strip().split();
                afmatrix.append([ float(snp[i]) for i in xrange(5,len(snp))]);
                alleles.append([snp[3],snp[4]]);
        File.close();
        return [afmatrix,alleles]

def read_cfmatrix(affile):
        afmatrix = []; alleles = [];
        File = open(affile,'r');
        for line in File:
                if line[0]=='#':continue;
                snp = line.strip().split();
                afmatrix.append([ float(snp[i]) for i in xrange(len(snp))]);
        File.close();
        return afmatrix


def compute_likelihood(aflist,genotypes,j,k):
        maxCLL = -1000000000; bestprop = 0;
        for t in range(0,101,5):
                prop1 = float(t)/100;
                CLL = 0;
                for i in range(len(aflist)):
                        genotype = genotypes[i];
                        maf1 = aflist[i][j]; maf2 = aflist[i][k];
                        p = maf1*prop1 + (1-prop1)*maf2;
                        if p < 0.00001: p = 0.00001;
                        elif p > 0.99999: p = 0.99999;
                        CLL += (genotype)*math.log(p) + (2.0-genotype)*math.log(1.0-p);
                        if genotype == 1: CLL += math.log(2);
                if CLL > maxCLL: maxCLL = CLL; bestprop = prop1;
                print prop1,CLL;
        print maxCLL,bestprop;


affile = sys.argv[1]; cffile = sys.argv[2]; pedfile = sys.argv[3];

[afmatrix,alleles] =read_matrix(affile);

cfmatrix = read_cfmatrix(cffile);
samples = len(cfmatrix); snps = len(afmatrix); pops = len(cfmatrix[0]);
print >>sys.stderr, 'samples:',samples, 'snps:',snps,'pops:',pops;

s=0;
File = open(pedfile,'r');
for line in File:
        G =line.strip().split();
        LL =0; genotypes = [];
        for j in xrange(snps):
                p1 =0.0;
                for k in xrange(pops):
                        p = afmatrix[j][k];
                        if p < 0.00001: p = 0.00001;
                        elif p > 0.99999: p = 0.99999;
                        p1 += cfmatrix[s][k]*p;

                if G[6+j*2] == alleles[j][0] and G[6+j*2+1] == alleles[j][0]: genotype =2;
                elif G[6+j*2] == alleles[j][0] and G[6+j*2+1] == alleles[j][1]: genotype =1;
                elif G[6+j*2] == alleles[j][1] and G[6+j*2+1] == alleles[j][0]: genotype =1;
                elif G[6+j*2] == alleles[j][1] and G[6+j*2+1] == alleles[j][1]: genotype =0;
                else: print >>sys.stderr, "Error genotype not defined";
                if genotype == 2: snpll = p1*p1;
                elif genotype ==1: snpll = 2*p1*(1.0-p1);
                elif genotype ==0: snpll = (1.0-p1)*(1.0-p1);
                LL += math.log(snpll)
                genotypes.append(genotype);
                #print 'SNP_NO',j,'af',afmatrix[j][0],'LL',LL,'sum',snpll;

        print G[0],G[1],cfmatrix[s],LL;

        compute_likelihood(afmatrix,genotypes,0,1);


        s +=1;
File.close();
