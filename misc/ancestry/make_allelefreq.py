#! /usr/bin/env python
# author Vikas Bansal vbansal@scripps.edu


import sys, os, glob, string,math, time
from subprocess import call

## read text file with allele frequency for each population, output allele-freq matrix 
def shared_allelefreqs(afile):
        File = open(afile);
        print >>sys.stderr, 'reading population allele frequencies';
        listofpopulations = {}; poptable = {}; snptable = {};
        lines = 0;
        for line in File:
                popsnp = line.strip().split();
                if popsnp[1] not in snptable: snptable[popsnp[1]] = [popsnp[2],int(popsnp[3]),popsnp[4],popsnp[7]];
                if popsnp[0] not in listofpopulations: listofpopulations[popsnp[0]] = 1;
                poptable[(popsnp[1],popsnp[0])] = popsnp[5];
                lines +=1;
        File.close();
        print >>sys.stderr, 'reading population allele frequencies';

        AFmatrix = []; snps =0;
        for snp in snptable.iterkeys():
                missing = 0; AF = [];
                for pop in listofpopulations.iterkeys():
                        try: AF.append(poptable[(snp,pop)]);
                        except KeyError: missing +=1;
                if missing ==0:
                        snpinfo = snptable[snp];
                        AFmatrix.append([snpinfo[0],snpinfo[1],snp,snpinfo[2],snpinfo[3],AF]);
                        snps +=1;


        AFmatrix.sort();
        print '#chrom','position','rsid','A1','A2',
        for pop in listofpopulations.iterkeys(): print pop,
        print;

        for i in xrange(snps):
                print AFmatrix[i][0],AFmatrix[i][1],AFmatrix[i][2],AFmatrix[i][3],AFmatrix[i][4],
                for f in AFmatrix[i][5]: print f,
                print;

        #print >>sys.stderr, 'noofpops',len(listofpopulations),'noofSNPs',lines/len(listofpopulations);
        return listofpopulations; # list of popID of all populations in the allele frequency file 


def convert_to_hg19(liftoverfile,afmatrixfile):
        File = open(liftoverfile);  LF = {};
        for line in File: snp= line.strip().split(); LF[snp[3]] = int(snp[1]);
        File.close();

        File = open(afmatrixfile);
        for line in File:
                if line[0] == '#': print line,
                else:
                        snp = line.strip().split();
                        try:
                                newposition = LF[snp[2]];
                                print snp[0],newposition,snp[2],
                                for i in xrange(3,len(snp)): print snp[i],
                                print;
                        except KeyError: pass;

        File.close();


# SNP_A-1780419   rs6576700       1       84647761        -       A       G
## make allele frequency file for Affy data 
def calculate_AF_strandfix(Affyfile,afile):

        RC = {}; RC['A'] ='T'; RC['T'] = 'A'; RC['G'] = 'C' ; RC['C'] = 'G';

        File = open(Affyfile); SNPS = {};
        for line in File:
                if line[0] == '#': continue;
                snp = line.strip().split();
                #print len(snp); 
                SNPS[snp[1]] = [snp[2],snp[3],snp[4],snp[5],snp[6]];
        File.close();

        File = open(afile);
        for line in File:
                snp = line.strip().split(',');
                if snp[0] == "rs_number": continue;
                try:
                        snpinfo = SNPS[snp[0]];
                        genotypes = snp[1:];
                        allele1 = snpinfo[3]; allele2 = snpinfo[4];
                        if snpinfo[2] == '-': allele1 = RC[allele1]; allele2 = RC[allele2];
                        counts = [0,0,0];
                        for i in xrange(len(genotypes)):
                                if genotypes[i][0] == allele1: counts[0] +=1;
                                elif genotypes[i][0] == allele2: counts[1] +=1;
                                else: counts[2] +=1;
                                if genotypes[i][1] == allele1: counts[0] +=1;
                                elif genotypes[i][1] == allele2: counts[1] +=1;
                                else: counts[2] +=1;
                        print snp[0],allele1,allele2,snpinfo,counts,genotypes;
                except KeyError: pass;


#shared_allelefreqs(sys.argv[1]);
calculate_AF_strandfix(sys.argv[1],sys.argv[2]);
