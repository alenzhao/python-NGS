#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler, random
from combinatorial import ncr, fet

#delta is used as index into list of pools, use delta = 0
def estimate_AF(allelecounts,pools,poolsize,SER,delta):  # poolsize is number of diploid individuals, not chromosomes 

	estimates = []; globalestimates = []; newestimates = []; totalalleles =0; psum = 0; wsum =0.000001;
	for p in xrange(pools): 
		if allelecounts[p+delta][2]+allelecounts[p+delta][3] < 10: continue;
		pi = float(allelecounts[p+delta][3]); # - SER*allelecounts[p+delta][2] -SER*allelecounts[p+delta][3]);
		pi /= allelecounts[p+delta][2]+allelecounts[p+delta][3]; #pi /= 1.000-2*SER;
		wi = float(2*(allelecounts[p+delta][2]+allelecounts[p+delta][3])); wi /= allelecounts[p+delta][2]+allelecounts[p+delta][3]+1;
		psum += pi*wi; wsum += wi; 
#		print pi,wi;
	pmean = psum/wsum; #print pmean;
	if pmean >= 0.9999: pmean = 0.999;
	if pmean <= 0.00001: pmean = 0.001;
		
	for p in xrange(pools):
		estimates.append([]); sum = -1000000; 
		for i in xrange(2*poolsize+1): 	
			af = float(i)/(2*poolsize);	P = af*(1-SER) + (1-af)*SER; 
			prior = ncr(2*poolsize,i) + math.log(pmean,10)*i + math.log(1-pmean,10)*(2*poolsize-i);
			#if p ==0: print prior;
			term = prior + math.log(P,10)*allelecounts[p+delta][3] + math.log(1-P,10)*allelecounts[p+delta][2]
			estimates[p].append(term);
			globalestimates.append(0.0); newestimates.append(0.0); totalalleles +=1;
			if term <= sum: sum += math.log(1+math.pow(10,term-sum),10); 
			else: sum = term + math.log(1+math.pow(10,sum-term),10);
			#print '%2d %2.1f ' %(i,estimates[p][i]),
		for i in xrange(2*poolsize+1): 	estimates[p][i] -= sum; 
#		print sum,
		#print '\n';

	for i in xrange(2*poolsize+1): globalestimates[i] = estimates[0][i]; 
	
	for p in xrange(1,pools):
		for r in xrange(0,(p+1)*(2*poolsize+1)): newestimates[r] = -10000 
		for q in xrange(0,2*poolsize+1):
			for r in xrange(0,p*(2*poolsize+1)):
				if newestimates[q+r] >= estimates[p][q] + globalestimates[r]:
					newestimates[q+r] += math.log(1+math.pow(10,estimates[p][q] + globalestimates[r]-newestimates[q+r]),10);
				else: 
					newestimates[q+r] = estimates[p][q] + globalestimates[r] + math.log(1+math.pow(10,newestimates[q+r]-estimates[p][q] - globalestimates[r]),10);
				#newestimates[q+r] += math.pow(10,estimates[p][q] + globalestimates[r]); 
		for r in xrange(0,(p+1)*(2*poolsize+1)):	globalestimates[r] = newestimates[r];
	
	min=0; lb =0; rb =0; bestll = -100000;
	for p in xrange(totalalleles): 
		if globalestimates[p] > bestll: min = p; bestll = globalestimates[p];
	lb = min; rb = min;
	while lb >= 0 and bestll-globalestimates[lb] <3: lb -=1;
	lb +=1;
	while rb <totalalleles and bestll-globalestimates[rb] <4: rb +=1;
	return [globalestimates,min,lb,rb];

			

def probabilisticFET(allelecounts,poolsH,poolsD,poolsize): # output from allele frequency estimation

	allelecounts.sort();
	HAF = estimate_AF(allelecounts,poolsH,poolsize/2,0.001,0);
	DAF = estimate_AF(allelecounts,poolsD,poolsize/2,0.001,poolsH);
	sum = -100000000;
	for i in xrange(HAF[2],HAF[3]):
		for j in xrange(DAF[2],DAF[3]):
			w = HAF[0][i] + DAF[0][j] + ncr(poolsize*poolsH,i) + ncr(poolsize*poolsD,j);
			if w <= sum: sum += math.log(1+math.pow(10,w-sum),10);
			else: sum = w + math.log(1+math.pow(10,sum-w),10);
	print sum,allelecounts;
	sum0 = sum;


	pvalue = 0;
	for perms in xrange(10000):
		random.shuffle(allelecounts);
		for i in xrange(poolsH): allelecounts[i][0] = 0; 
		for i in xrange(poolsD): allelecounts[i+poolsH][0] = 1; 
		allelecounts.sort();
		HAF = estimate_AF(allelecounts,poolsH,poolsize/2,0.001,0);
		DAF = estimate_AF(allelecounts,poolsD,poolsize/2,0.001,poolsH);
		sum = -100000000;
		for i in xrange(HAF[2],HAF[3]):
			for j in xrange(DAF[2],DAF[3]):
				w = HAF[0][i] + DAF[0][j] + ncr(poolsize*poolsH,i) + ncr(poolsize*poolsD,j);
				if w <= sum: sum += math.log(1+math.pow(10,w-sum),10);
				else: sum = w + math.log(1+math.pow(10,sum-w),10);
		if sum < sum0: pvalue +=1;
		if pvalue >= 2: break;
		print sum;

	return float(pvalue+1)/(perms+2);
	print 'pvalue',float(pvalue+1)/(perms+2);
	
# this uses the old VCF format of CRISP, now changed 
def calculate_FETpvalues(vcffile,poolDX,poolsize):
	File = open(vcffile);
	for line in File:
		if line[0] == '#' and line[1] == '#': continue;
		variant = line.strip().split();
		if variant[0] == '#CHROM': 
			samples = len(variant)-9; 
			#samplelist = variant[9+offset:]; samples = len(samplelist);
			#print samplelist,samples;
			continue;
		
		if 'SNP' in variant[2] or 'INDEL' in variant[2]: 
			H0 = 0; H1 = float(0); D0 = 0; D1 = float(0);
			allelecounts = []; allelecountsD = []; poolsH = 0; poolsD = 0;
			for i in xrange(samples):
				if poolDX[i] == -1: continue;  # ignore pool for case control analysis 

				counts = variant[i+9].split(':');
				total = int(counts[0]) + int(counts[1]);
				#if total < 240: continue;
				alt = float(counts[1])*poolsize; alt /= (total+0.01); 
				#if alt < 0.5: alt = 0;  
				#if alt > 0.5 and alt < 1: alt = 1;
				
				if poolDX[i] == 0: 
					H0 += poolsize; H1 += alt; 
					allelecounts.append([0,0,int(counts[0]),int(counts[1])]); poolsH +=1;
					# 0/1 bit of allelecounts stores case-control status...
				elif poolDX[i] ==1: 
					D0 += poolsize; D1 += alt; 
					allelecounts.append([1,1,int(counts[0]),int(counts[1])]); poolsD +=1;
			#	print counts,total,alt;
			if H1 + D1 >= 5:
				#print allelecounts;
				allelecounts.sort();
				HAF = estimate_AF(allelecounts,poolsH,poolsize/2,0.001,0);
				DAF = estimate_AF(allelecounts,poolsD,poolsize/2,0.001,poolsH);
				"""
				for p in xrange(HAF[2],HAF[3]): print '%3d %2.2f ' %(p,HAF[0][p]),
				print 'Binomial maxll',HAF[1],HAF[0][HAF[1]];
				for p in xrange(DAF[2],DAF[3]): print '%3d %2.2f ' %(p,DAF[0][p]),
				print 'Binomial maxll',DAF[1],DAF[0][DAF[1]];
				"""

				pvalue = fet(int(round(H1,0)),int(H0),int(round(D1,0)),int(D0));
				if pvalue[0] < -2: 
					pvalueperm = probabilisticFET(allelecounts,poolsH,poolsD,poolsize); # output from allele frequency estimation
					print 'PERM',
				else: 
					pvalueperm = 1; 		
					print 'FET',
				print math.log(pvalueperm,10),pvalue[0],variant[0],variant[1],variant[2],variant[3],variant[4],variant[5],variant[6],variant[7],
				print '%2.1f %2.1f %2.1f %2.1f' %(H0,H1,D0,D1),#'Control',float(H1)/(H0+0.001),'Case',float(D1)/(D0+0.001),
				if pvalue < 0.001: print 'LOW';
				else: print;
			


