#! /usr/bin/env python
# AUTHOR VIKAS BANSAL last edited sept 7 2012, vbansal@scripps.edu 

# [vbansal-scripps@tcc-5-19 VariantCalls-aug252012]$ python ../Association/pythoncode/genetests.py CRISP.genes/variants.SRR.vcf ../Phenotypes-pools/192pools.phenotypes 40

### code for doing gene based association tests for T2D sequencing data ##

## python ~/PROGRAMS/JOINTCODE-coral/PYTHON-scripts/NGS-library-python/genetests-triallelic.py CRISP.genes.VCF.header.sorted ../Phenotypes-pools/192pools.phenotypes 40 > triallelic.pvalues

import sys, os, glob, string, subprocess,time, math, re, compiler, random
#from scipy.stats import chi2  # chi square test p-value 
from combinatorial import ncr, fet
from associationfunctions import estimate_AF,probabilisticFET,calculate_FETpvalues ## old functions in separate file 

def read_phenotypefile(statusfile):
	## read the case control status file  ### 
	poolDX = {}; casepools = 0; controlpools = 0; totalpools =0;
	File = open(statusfile);
	for line in File: 
		pool = line.strip().split(); 
		poolDX[pool[0]] = int(pool[1]);
		if int(pool[1]) == 0: controlpools +=1; 
		if int(pool[1]) == 1: casepools +=1; 
		totalpools +=1; 
	File.close();
	print >>sys.stderr, "totalpools",totalpools,"CASEpools",casepools,"controlpools",controlpools;
	return poolDX; 

############################################################################################
# this uses the new VCF format of CRISP
def calculate_FETpvalues_new(vcffile,poolstatustable,poolsize):
	maxMAF = 0.01; # maximum minor allele frequency, variants above this are not used
	poolDX = [];  # phenotype status for each pool
	Variants = []; 

	File = open(vcffile);
	for line in File:
		if line[0] == '#' and line[1] == '#': continue;
		variant = line.strip().split('\t');
		if variant[0] == '#CHROM': 
			for i in xrange(9,len(variant)):
				sampleid = variant[i].split('/')[-1].split('.')[0];
				#print sampleid;
				try: status = poolstatustable[sampleid]; 
				except KeyError: status= -1; 
				poolDX.append(status);
				print >>sys.stderr, status,
			#samples = len(variant)-9; print 'samples',samples,len(variant);
			print >>sys.stderr;
			continue;

		samples = len(variant)-9;
		chrom = variant[0]; position = int(variant[1]); refallele = variant[3]; varalleles = variant[4].split(','); 

		if len(varalleles) ==2: triallelic=1;
		else: triallelic=0;
		if len(varalleles) > 2: continue; # ignore multi-allelic variants for now

		#if triallelic ==0 or position != 64527465: continue;

		H0 = 0.0; H1 = 0.0; H2=0.0; D0 = 0.0; D1 = 0.0; D2=0.0;
		casepools=0; controlpools=0;
		for i in xrange(samples):
			if poolDX[i] == -1: continue;  # ignore pool for case control analysis 
			try: 
				counts = variant[i+9].split(':');
				readsf = counts[2].split(','); readsr = counts[3].split(','); 
				total = int(readsf[0]) + int(readsf[1]) + int(readsr[0]) + int(readsr[1]); 
				if triallelic==1: total += int(readsf[2])+int(readsr[2]); # tri-allelic

				alt = (float(readsf[1]) + float(readsr[1]))*poolsize; alt /= (total+0.0001); 
				if triallelic==1: alt2 = (float(readsf[2]) + float(readsr[2]))*poolsize; alt2 /= (total+0.0001); 
				else: alt2=0.0;

				sf = float(2*total)/(2*total + poolsize); 
				if poolDX[i] == 0: 
					#H0 += poolsize; H1 = alt; 
					H0 += sf*float(poolsize); H1 += sf*alt; H2 += sf*alt2; controlpools +=1;
					# 0/1 bit of allelecounts stores case-control status...
				elif poolDX[i] >= 1: 
					D0 += sf*float(poolsize); D1 += sf*alt; D2 += sf*alt2; casepools +=1; 
					#D0 += poolsize; D1 += alt; 
			except IndexError: print 'Exception',i,samples,variant; 

		#if H1 + D1 < 10 or (H0-H1 + D0-D1 <10): lowcountvariants +=1; continue; 
		pvalue = [0,0,0]; pvalue2=[0,0,0];
		if H1 + D1 >= 4: pvalue = fet(int(round(H1,0)),int(round(H0,0)),int(round(D1,0)),int(round(D0,0)));
		print '%0.2f %.2f:%.2f %.2f:%.2f'  %(pvalue[0],H0,H1,D0,D1),#'Control',float(H1)/(H0+0.001),'Case',float(D1)/(D0+0.001),
		if triallelic ==1 and H2+D2 >=4: 
			pvalue2 = fet(int(round(H2,0)),int(round(H0,0)),int(round(D2,0)),int(round(D0,0)));
			print '2nd-allele %0.2f %.2f:%.2f %.2f:%.2f'  %(pvalue2[0],H0,H2,D0,D2),

		print '%3s %9s %10s %10s %20s ' %(variant[0],variant[1],variant[3],variant[4],variant[5]+':'+variant[6]),
		print variant[7];

############################################################################################

if len(sys.argv) < 4: 
	print >>sys.stderr, "python pooled-association.py VCFfile poolDXstatusfile poolsize(haploid)"; sys.exit();

random.seed();
vcffile=sys.argv[1];  statusfile = sys.argv[2]; poolsize = int(sys.argv[3]);

poolDX = read_phenotypefile(statusfile); 
calculate_FETpvalues_new(vcffile,poolDX,poolsize);



############################################################################################

	
	



