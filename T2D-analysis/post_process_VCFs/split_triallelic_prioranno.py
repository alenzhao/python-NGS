#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler

#last modified jan 23 2014
## python filter_rarevariants_VCF.py 6pops.allelefreq 190pools.ancestry.5pops combined.vcf.annotated > b

## CT pvalues can also be split, EMfail/EMpass...

if len(sys.argv) < 2: print >>sys.stderr, "python split_triallelic.py combined.vcf.annotated 0/1/2(printVCF)"; sys.exit();

printVCF=0;
if len(sys.argv) > 2: printVCF = int(sys.argv[2])

variants=0; overlapping=0; poollist = [];
File = open(sys.argv[1]); # pooled VCF file 
for line in File:
	if line[0] == '#' and line[1] == 'C' and line[2] == 'H': var = line.strip().split('\t'); poollist = var;

	if line[0] == '#': 
		if printVCF ==1: print line,
		continue;
	var = line.strip().split('\t'); chrom = var[0]; position = int(var[1]); alleles = var[4].split(','); 

	splitvariant =0; counts = [0,0]; pools = 0;

	if len(alleles) ==2: ## only split tri-allelic, variants with 4 alleles are even more difficult..
		info = var[7].split(';'); pools = len(var)-9;
		for i in xrange(len(info)): 
			if info[i].split('=')[0] == 'AC': counts = [int(info[i].split('=')[1].split(',')[0]),int(info[i].split('=')[1].split(',')[1])];

		if (counts[0]*20 <= pools and counts[1]*20 <= pools) or counts[0] <= 5 or counts[1] <= 5: splitvariant = 1; 

	if splitvariant ==0:
		if printVCF==1: print line, 
		elif printVCF==2:
			var[2] = '.';
			for i in xrange(8): sys.stdout.write(var[i] + "\t")
			sys.stdout.write('\n')
			
		else: pass
	else:
		if printVCF ==0: print counts,alleles[0],var[3],alleles[1],var[0:8];
		if printVCF ==2: var[2] = '.';
		if printVCF >= 1:
			for i in xrange(4): sys.stdout.write(var[i] + "\t")
			sys.stdout.write(alleles[0] + '\t')
			sys.stdout.write(var[5].split(',')[0] + '\t')
			for i in xrange(6,7): sys.stdout.write(var[i] + '\t')
			sys.stdout.write(var[7] + ";FLAG=SPLIT\t");

			if printVCF != 2: 
				sys.stdout.write(var[8]);  # needs to be split as well...
				for i in xrange(9,len(var)):
					genotype = var[i].split(':'); 
					newgenotype = genotype[0].split(',')[0]; 
					newvar = ':'.join([newgenotype] + genotype[1:]);
					sys.stdout.write('\t' + newvar);
					#if int(newgenotype) > 0: print "\t%s" %(newvar),

			sys.stdout.write('\n')

		if printVCF >= 1:
			for i in xrange(4): sys.stdout.write(var[i] + "\t")
			sys.stdout.write(alleles[1]+ '\t')
			sys.stdout.write(var[5].split(',')[1] + '\t')
			for i in xrange(6,7): sys.stdout.write(var[i] + '\t')
			sys.stdout.write(var[7] + ";FLAG=SPLIT\t");

			if printVCF != 2: 
				sys.stdout.write(var[8]);
				for i in xrange(9,len(var)):
					genotype = var[i].split(':'); 
					newgenotype = genotype[0].split(',')[1]; 
					newvar = ':'.join([newgenotype] + genotype[1:]);
					sys.stdout.write('\t' + newvar);
			sys.stdout.write('\n')
	
		if printVCF ==0: print counts[0],counts[1];
		
		#print var[0],var[1],var[3],alleles,var[2],var[7];


File.close();
