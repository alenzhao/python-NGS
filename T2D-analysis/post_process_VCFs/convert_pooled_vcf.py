#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler

## coded july 15 2015, code to convert CRISP pooled allele count to genotype 0/0/0/0/0/1 format so that it is compatible with standard tools

if len(sys.argv) < 2: print >>sys.stderr, "python split_triallelic.py combined.vcf.annotated 0/1/2(printVCF)"; sys.exit();

variants=0; 
poolsize = -1; 
if len(sys.argv) > 2: poolsize = int(sys.argv[2]); 

File = open(sys.argv[1]); # pooled VCF file 

for line in File:

	if line[0] == '#': 
		print line,
		continue;

	var = line.strip().split('\t'); chrom = var[0]; position = int(var[1]); alleles = var[4].split(',');

	"""
	"""
	for i in xrange(7): sys.stdout.write(var[i] + "\t")
	sys.stdout.write(var[7] + ";FLAG=SPLIT\t");
	#sys.stdout.write(var[8]);
	sys.stdout.write('GT:GQ:ADf:ADr:ADb'); # replace AlleleCount by Genotype 

	for i in xrange(9,len(var)):
		genotype = var[i].split(':'); counts = genotype[0].split(','); 

		if poolsize > 0: 
			refsum = poolsize;
			for c in xrange(len(counts)): refsum -= int(counts[c]); 
			GV = [];
			for c in xrange(refsum): GV.append('0'); 
			for a in xrange(len(counts)): 
				for c in xrange(int(counts[a])):  GV.append(`(a+1)`); 
		#print 'genotype','/'.join(GV);					
		else: ## variable pool size so genotype is different 
			GV = [];
			for a in xrange(len(counts)): 
				for c in xrange(int(counts[a])):  GV.append(`(a)`); 
			
			

		sys.stdout.write('\t' + '/'.join(GV) + ':' +  ':'.join(genotype[1:]));
	
	sys.stdout.write('\n')


File.close();
