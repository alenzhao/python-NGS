#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler
space = re.compile(r'\s+');
compiler.parseFile(sys.argv[0]); 

mincov = 5;

def extract_snps(snipseqfile):
	currlocus = ""; currpos = -10; genotypelist = [];
	File=open(snipseqfile,'r');  reads =0;
	for s in File:
		line = space.split(s);
		if line[0] != "GT" or line[5] != "Y": continue;
		genotype = int(line[9]);allele1 = line[13]; allele2 = line[14];
		scores = [float(line[35]),float(line[36]),float(line[37])];
		if genotype ==0: G = allele1+allele1; conscore = int((scores[0]-scores[1])*10);
		elif genotype ==1: G = allele1 + allele2; conscore = min( int((scores[1]-scores[0])*10), int((scores[1]-scores[2])*10));
		elif genotype ==2: G = allele2 + allele2; conscore = int((scores[2]-scores[1])*10); 
		else: G = 'NN'; conscore = 0; 
		coverage = int(line[24]);
		if coverage < mincov: G = 'NN'; conscore =0;

		print '%9s %10s %1s %15s %2s %4d |  %3d ALLELES %2d %2d %2d %2d %2d %2d %2d %2d Likelihoods %+2.2f %+2.2f %+2.2f' %(line[1],line[2],allele1,line[3],G,conscore,coverage,int(line[16]),int(line[17]),int(line[18]),int(line[19]),int(line[20]),int(line[21]),int(line[22]),int(line[23]),scores[0],scores[1],scores[2]),
		if line[11] == '1': print 'dbSNP';
		else: print 'novel';

	File.close();

extract_snps(sys.argv[1]);


