#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler
space = re.compile(r'\s+');
compiler.parseFile(sys.argv[0]);

# python pipeline_coverage.py 080703_HWI-EAS76_0004/MAQ_SNP_2008-09-09-193705/5.J400249.Adapter1.pileup genomes/Sanofi.bed 
BEDFILE =1;

if len(sys.argv) < 3:
        print "please enter a pileup file and .bed file ";
        sys.exit();
else:
        pileupfile = sys.argv[1];
        ampliconfile = sys.argv[2]; #.bed file from genomes/

#############################READ .pileup file for storing coverage stats          ###############################

index_table = {}; File = open(pileupfile,'r');
for line in File: A = line.split(); index_table[(A[0],int(A[1]))] = int(A[3]); prev = int(A[1]);
File.close();

######################### now compute statistics on coverage per amplicon and output to file #####################

File = open(ampliconfile,'r');
for s in File:
        line = space.split(s);
        chrom = line[0]; start = int(line[1]); end = int(line[2]);
        totalcov =0; positions =0; cov =0; covsq = 0; mean =0; median =0; std=0;
        for i in range(start,end+1):
                try:
                        cov = index_table[(chrom,i)];
                        totalcov += cov; covsq += cov*cov;
                except KeyError: pass;

        positions = end-start+1;
        mean = float(totalcov)/float(positions); std = float(covsq)/float(positions) - mean*mean;
        print '%6s %12s %12s\t%5d\t%5d' % (line[0],line[1],line[2],int(mean),int(math.sqrt(std)))
File.close();
