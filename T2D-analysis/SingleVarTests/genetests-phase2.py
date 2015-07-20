#! /usr/bin/env python
# AUTHOR VIKAS BANSAL 
### code for doing gene based association tests for T2D sequencing data ## feb 6 2015

## edited to work for phase2, no DiabeticComplication data, CRISP output is slightly different (no meanAC) 

## exclude pools with low coverage per site for calculating p-value 

## python genetests-phase2.py ../phase2-132pools/132pools.CRISP.082413.variants/132pools.CRISP.082413.combined.annotated.vcf ../phase2-132pools/132pools.phenotypes 48 > pvalues.phase2


## outputs info on early onset allele frequencies etc...

import sys, os, glob, string, subprocess,time, math, re, compiler, random
#from scipy.stats import chi2  # chi square test p-value 
from combinatorial import ncr, fet
from associationfunctions import estimate_AF,probabilisticFET,calculate_FETpvalues ## old functions in separate file 

PRINTGENOTYPES =0;
USE_MEAN_AC = 0;

## read the phenotype file with poolID (G1P54), 2,3 etc
def read_phenotypefile(statusfile):
	## read the case control status file  ### 
	poolDX = {}; casepools = 0; controlpools = 0; totalpools =0; earlyonset=0; DC=0;
	File = open(statusfile);
	for line in File: 
		## allow multiple phenotype codes for each pool separated by commalist 
		pool = line.strip().split(); phenotype = pool[1].split(',');  
		if len(phenotype) ==1: poolDX[pool[0]] = [int(phenotype[0])]; 
		elif len(phenotype) ==2: poolDX[pool[0]] = [int(phenotype[0]),int(phenotype[1])]; 
		elif len(phenotype) ==3: poolDX[pool[0]] = [int(phenotype[1]),int(phenotype[2]),int(phenotype[3])]; 

		if int(phenotype[0]) == 0: controlpools +=1; 
		if int(phenotype[0]) >= 1: casepools +=1; 
		if int(phenotype[0]) == 2: earlyonset +=1; 
		if len(phenotype) >=2 and int(phenotype[1]) == 3: DC +=1; 
		totalpools +=1; 
	File.close();
	print >>sys.stderr, "totalpools",totalpools,"CASEpools",casepools,"controlpools",controlpools,'earlyonset',earlyonset,'T2D-complication',DC;
	return poolDX; 


def determine_poolphenotype(vcffile,poolstatustable):
	poolDX = [];  # phenotype status for each pool
	File = open(vcffile);
	for line in File:
		if line[0] == '#' and line[1] == '#': continue;
		variant = line.strip().split('\t');
		if variant[0] == '#CHROM': 
			for i in xrange(9,len(variant)):
				sampleid = variant[i].split('/')[-1].split('.')[0];
				#print sampleid;
				try: status = poolstatustable[sampleid]; 
				except KeyError: status= [-1]; 
				poolDX.append(status);
				print >>sys.stderr, status,
			#samples = len(variant)-9; print 'samples',samples,len(variant);
			print >>sys.stderr;
			continue;
	File.close();
	return poolDX; 


############################################################################################
# this uses the VCF format where we have allele counts (ML, MEAN, VARIANCE) instead of readcounts
# ignore multi-allelic variants for now
def calculate_FETpvalues_allelecounts(vcffile,poolstatustable,poolsize):
	poolDX = determine_poolphenotype(vcffile,poolstatustable);  # phenotype status for each pool
	variants =0; trivariants =0;

	File = open(vcffile);
	for line in File:
		if line[0] == '#': continue;
		variant = line.strip().split('\t');
		samples = len(variant)-9;
		chrom = variant[0]; position = int(variant[1]); refallele = variant[3]; varalleles = variant[4].split(','); 

		if len(varalleles) ==2: triallelic=1; trivariants +=1; 
		else: triallelic=0; variants +=1; 

		if len(varalleles) >= 3 or len(varalleles) ==4: 
			#print >>sys.stderr, '##triallelic',chrom,position,refallele,varalleles,variant[2];
			continue; # ignore multi-allelic variants for now

		## healthy (controls), D (disease), E (early onset), DC (diabetic complications)  
		H0 = 0.0; H1 = 0.0; H2 = 0.0; D0 = 0.0; D1 = 0.0; D2 = 0.0; E0 = 0.0; E1= 0.0; E2 = 0; DC0 = 0; DC1 = 0;  DC2 = 0;
		casepools=0; controlpools=0;
		for i in xrange(samples):
			if poolDX[i][0] == -1: continue;  # ignore pool for case control analysis 
			try: 
				genotypes = variant[i+9].split(':');
				if triallelic==0: 
					MLAC = int(genotypes[0]); #meanAC = float(genotypes[2]); QAC = float(genotypes[1]); 
					MLAC2 = 0;
				else: 
					MLAC = int(genotypes[0].split(',')[0]); #meanAC = float(genotypes[2]); QAC = float(genotypes[1]); 
					MLAC2 = int(genotypes[0].split(',')[1]);

				#if genotypes[3] == '-inf': varAF = -0.1*QAC; 
				#else: varAF = float(genotypes[3]);
 
				if poolDX[i][0] == 0: H0 += poolsize; H1 += MLAC; H2 += MLAC2;
				elif poolDX[i][0] >= 1:  D0 +=poolsize; D1 += MLAC; D2 += MLAC2;
				if poolDX[i][0] == 2: E0 += poolsize; E1 += MLAC; E2 += MLAC2;
				if len(poolDX[i]) >=2 and poolDX[i][1] == 3: DC0 += poolsize; DC1 += MLAC; DC2 += MLAC2; 

			except IndexError: print 'Exception',i,samples,variant; 

		## calculate p-value between cases and controls 
		if (H1 + D1 >= 4): pvalue = fet(int(round(H1,0)),int(round(H0,0)),int(round(D1,0)),int(round(D0,0)));
		else: pvalue = [0,0,0];

		## calculate p-values between early onset and controls
		if (H1+E1 >=4): pvalue1 = fet(int(round(H1,0)),int(round(H0,0)),int(round(E1,0)),int(round(E0,0)));
		else: pvalue1 = [0,0,0]; 

		## calculate p-value between early onset and late-onset
		if (D1+E1 >=4): pvalue2 = fet(int(round(D1,0)),int(round(D0,0)),int(round(E1,0)),int(round(E0,0)));
		else: pvalue2 = [0,0,0]; 
		

		print '%0.2f %0.2f %.1f/%.1f %.1f/%.1f %.1f/%.1f'  %(pvalue[0],pvalue1[0],H1,H0,D1,D0,E1,E0),#'Control',float(H1)/(H0+0.001),'Case',float(D1)/(D0+0.001),
		#print '%0.2f %0.2f' %(pvalue1[0],pvalue2[0]),
		#print '%0.4f %0.4f %0.4f %0.4f' %(H1/H0,D1/D0,E1/E0,DC1/DC0),
		print '%0.4f %0.4f %0.4f' %(H1/H0,D1/D0,E1/E0),

		if triallelic ==1: 
			if H2 + D2 >=4: pvalue_tri = fet(int(round(H2,0)),int(round(H0,0)),int(round(D2,0)),int(round(D0,0)));
			else: pvalue_tri = [0,0,0];
			print 'TRIALLELIC:%0.2f:%.1f:%.1f'  %(pvalue_tri[0],H2,D2),
		else: print '-',

		print '%3s %9s %s %10s %10s %8s %10s' %(variant[0],variant[1],variant[2],variant[3],variant[4],variant[5],variant[6]),
		print variant[7], 

		if PRINTGENOTYPES ==1: 
			for i in xrange(samples): 
				if poolDX[i][0] == -1: continue;  # ignore pool for case control analysis 
				genotypes = variant[i+9].split(':');
				print variant[i+9],
				#MLAC = int(genotypes[0]); meanAC = float(genotypes[2]); QAC = float(genotypes[1]); 
				#if genotypes[3] == '-inf': varAF = -0.1*QAC; 
				#else: varAF = float(genotypes[3]); 
				#print '%2d:%0.2f:%0.2f' %(MLAC,meanAC,math.sqrt(pow(10,varAF))),
		print;
	print >>sys.stderr, "variants evaluated",variants,"triallelic or more",trivariants;


############################################################################################
TYPE = 0;

if len(sys.argv) < 4: 
	print >>sys.stderr, "python pooled-association.py VCFfile phenotypefile poolsize allelecounts/readcounts"; sys.exit();
if len(sys.argv) >5: PRINTGENOTYPES = int(sys.argv[5]); print >>sys.stderr, "genotypes will be printed";
if len(sys.argv) >4: TYPE = int(sys.argv[4]);

vcffile=sys.argv[1];  statusfile = sys.argv[2]; poolsize = int(sys.argv[3]);
poolDX = read_phenotypefile(statusfile); 

if TYPE ==0: calculate_FETpvalues_allelecounts(vcffile,poolDX,poolsize);
else: calculate_FETpvalues_readcounts(vcffile,poolDX,poolsize);



############################################################################################

	
	



