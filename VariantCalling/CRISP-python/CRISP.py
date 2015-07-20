#! /usr/bin/env python

###############################################################################################################################
# CRISP: Comprehensive Read Analysis for Identification of SNPs (and short indels) from Pooled sequencing data
# Author: Vikas Bansal, vbansal@scripps.edu 
# (c) 2010 Vikas Bansal, Scripps Genomic Medicine 

# This program is distributed free for academic/research use only. 
# This program comes with no guarantees. Use it at your own risk. 
# This program has been tested using python version 2.4 on a linux x86_64 system and may not work on other platforms. 
# This program may not be redistributed in any form without permission of the author. 
# Permission is granted to modify it for personal use but all modified versions should include this notice.

###############################################################################################################################


import sys, os, glob, string, subprocess,time, math, random 
from optparse import OptionParser
from array import *
from scipy import stats


def ncr(n,r):
	ll = 0;
	for i in xrange(min(r,n-r)): ll += math.log(float(n-i)/(i+1),10);
	return ll;

def mcr(n,rlist):
	ll =0; k=n;
	for i in xrange(len(rlist)):
		for j in xrange(rlist[i]): ll += math.log(k,10)-math.log(j+1,10); k -=1;
	return ll;

def fet(c1,r1,c2,r2):
	minc =0;
	if r2 < c1+c2: minc = c1+c2-r2;
	maxc = c1+c2;
	#if r1 < c1+c2: maxc = c1+c2-r1;
	if r1 < c1+c2: maxc = r1;
	ptable = ncr(r1,c1) + ncr(r2,c2); 	

	c11 = minc; c21 = c1+c2-minc;
	p1 = ncr(r1,c11); p2 =  ncr(r2,c21); 
	pvalue0 = 0; pvalue1 =0; sum =0; 

	max = 0; 
	while c11 <= maxc:
		if p1 + p2 > max: max = p1+p2;
		if r1-c11 > 0: p1 += math.log(r1-c11,10) - math.log(c11+1,10);
		if c21 > 0: p2 += math.log(c21,10) - math.log(r2-c21+1,10);
		c11 +=1; c21 -= 1;

	c11 = minc; c21 = c1+c2-minc;
	p1 = ncr(r1,c11); p2 =  ncr(r2,c21); 
	while c11 <= maxc:
		val = math.pow(10,p1+p2-max); sum += val;
		if p1+p2 <= ptable: pvalue0 += val;
		#print p1+p2,val,pvalue0,sum; print 'TABLE',c11,r1,c21,r2,p1+p2,ptable;
		if r1-c11 > 0: p1 += math.log(r1-c11,10) - math.log(c11+1,10);
		if c21 > 0: p2 += math.log(c21,10) - math.log(r2-c21+1,10);
		c11 +=1; c21 -= 1;
		#p1 = ncr(r1,c11); p2 =  ncr(r2,c21); 
	if sum > 0: return pvalue0/sum;
	else: return pvalue0/(sum+0.01); 
			

def compute_likelihood(ctable):
	ll = 0;
	for i in range(len(ctable)): 	
		if len(ctable[i]) ==2: r = ctable[i][1]; n = ctable[i][0];
		if len(ctable[i]) ==4: r = ctable[i][1]+ctable[i][3]; n = ctable[i][0]+ctable[i][2];
		if r ==0 or n ==r: continue;
		elif r ==1 or r == n-1: ll += math.log(n,10);
		elif r ==2 or r == n-2: ll+= math.log(n*(n-1)*0.5,10);
		else: 
			ll0 = 0;
			for i in range(min(r,n-r)): ll0 += math.log(float(n-i)/(i+1),10);
			ll += ll0;
	return ll;


def correctionfactor(ctable):
	sizes = {};
	for i in xrange(len(ctable)): 
		try: sizes[(ctable[i][0],ctable[i][1])] +=1;
		except KeyError: sizes[(ctable[i][0],ctable[i][1])] =1;
	slist = [];
	for b in sizes.iterkeys(): slist.append([b[0],b[1],sizes[b]]);
	slist.sort();

	#print 'slist',slist;
	factor = 0;
	i=0; s= slist[i][0]; total = 0; rlist = [];
	while i < len(slist):
		if slist[i][0] == s: 	rlist.append(slist[i][2]); total += slist[i][2]; 
		else:
			factor += mcr(total,rlist);
			s = slist[i][0]; total = slist[i][2]; rlist = [slist[i][2]];
		i +=1;
	#print rlist,total;
	factor += mcr(total,rlist);
	#print 'CORR',slist, factor;
	return factor;


def table_constants(ctable):
	n = len(ctable);
	ones = 0; reads = 0;
	for i in xrange(n): ones += ctable[i][1]; reads += ctable[i][0];
	table_index = array('I'); offset = 0;
	for i in xrange(n): 
		for j in range(ctable[i][0]): table_index.append(i); 
		offset += ctable[i][0];
	return [reads,ones,table_index];


def simulate_random_table_efficient(ctable,reads,ones,table_index,new_ctable): # new_ctable is already there 
	selected_list = {}; newll = 0;
	for i in xrange(len(ctable)): new_ctable[i][1] = 0; 
	for i in xrange(ones):
		r = int(random.random()*reads);
		while r in selected_list: r = int(random.random()*reads);
		selected_list[r] = 1; b = table_index[r]; 
		newll += math.log(float(new_ctable[b][0]-new_ctable[b][1])/(new_ctable[b][1]+1),10); new_ctable[b][1] +=1; 
	return newll;


def pvalue_chi2(table):
	# table is 2 x k matrix
	reads =0; ones =0;
	for i in range(len(table)): reads += table[i][0]; ones += table[i][1];
	statistic =0;
	for i in xrange(len(table)):
		Ei0 = table[i][0]*(reads-ones)/reads;
		Ei1 = table[i][0]*ones/reads;
		if Ei0 == 0 or Ei1 ==0: return 1;
		d = 0.0;
		statistic += (abs(table[i][0]-table[i][1]-Ei0)-d)*(abs(table[i][0]-table[i][1]-Ei0)-d)/Ei0;
		statistic += (abs(table[i][1]-Ei1)-d)*(abs(table[i][1]-Ei1)-d)/Ei1;
	degreesFreedom = len(table)-1;
	return stats.chisqprob(statistic,degreesFreedom)
	return 1;


def compute_pvalue_exact(ctable,R,A,prob,NCRtable):
	if A > R: return -1;
	if len(ctable) ==1: 
		if ctable[0][0] > 0: cp = NCRtable[(ctable[0][0],A)];
		else: cp = 0;
		if cp > prob:  return -1;
		return cp;
	else:
		sum =0;
		for i in range(min(ctable[0][0]+1,A+1)):  # correct bounds for i  
			if ctable[0][0] > 0: cp =NCRtable[(ctable[0][0],i)];
			else: cp = 0; 
			a = compute_pvalue_exact(ctable[1:],R-ctable[0][0],A-i,prob-cp,NCRtable);
			if a > -1: sum += math.pow(10,a + cp-prob);
		if sum > 0: return math.log(sum,10) + prob;
		return -1; 
	
def compute_NCRtable(maxn):
	NCRtable = {};
	for n in range(1,maxn+1):
		value = 0; NCRtable[(n,0)] = 0; #print n,0,0;
		for r in range(1,n): value += math.log(float(n-r+1)/(r),10); NCRtable[(n,r)] = value; # print n,r,value;
		NCRtable[(n,n)] = 0; # print n,n,0;
	return NCRtable;

def compute_pvalue_permutation(a,cl,stranded,iter):
	pvalue = 0; pvalue2 = 0;
	newtable = []; 
	for k in xrange(len(a)): newtable.append([a[k][0],0]); 
	[reads,ones,table_index] = table_constants(a);
	ifactor = correctionfactor(a);

	for k in xrange(1,iter+1): 
		for i in xrange(len(a)): newtable[i][1] = 0; 
		clra =0;
		for i in xrange(ones):
			j = int(random.random()*(reads-i)) + i; b = table_index[j]; table_index[j] = table_index[i]; table_index[i] = b;
			clra += math.log(float(newtable[b][0]-newtable[b][1])/(newtable[b][1]+1),10); newtable[b][1] +=1; 
	#	clra = simulate_random_table_efficient(a,reads,ones,table_index,newtable); 
		factor = correctionfactor(newtable);
#		print newtable,'|',cl,ifactor,clra,factor,pvalue,'|';
		if clra + factor <= cl+ifactor: 	pvalue +=1; 
		if k <= 20 and pvalue >= 3: break;
		if k <= 100 and pvalue >= 5: break;
		if k <= 1000 and pvalue >= 5: break;
	#print pvalue,k;
	if pvalue ==0: return float(pvalue+1)/(k+1);
	else: return float(pvalue)/k;

def compute_pvalue_fisher(a,cl,stranded,iter,NCRtable):
	[reads,ones,table_index] = table_constants(a);
	maxcol = 0;
	for i in range(len(a)): 
		if a[i][0] > maxcol: maxcol = a[i][0];
	if ones < 70 and maxcol < 1000 and len(a) > 2 and len(a) < 8:
		p = compute_pvalue_exact(a,reads,ones,cl,NCRtable);
		#print 'final',cl,ncr(reads,ones),reads,ones,p,p-ncr(reads,ones);  
		return math.pow(10,p-ncr(reads,ones));
	else:
		p= compute_pvalue_permutation(a,cl,stranded,iter);
		return p;

def minreads_pvalue(R,A,e):
	e1 = math.log(e,10); e2 = math.log(1-e,10); ll = R*e2; pvlog = ll; sum = ll; 
	for r in range(1,A+1): 
		pvlog += math.log(R-r+1,10) - math.log(r,10) + e1-e2;
		sum += math.log(1+math.pow(10,pvlog-sum),10);
	return sum;

def binomial_pvalue(R,A,e):
	e1 = math.log(e,10); e2 = math.log(1-e,10); 
	ll = ncr(R,A) + A*e1 + (R-A)*e2; pvlog = ll; sum = ll;
	for r in range(A+1,R+1): 
		pvlog += math.log(R-r+1,10) - math.log(r,10) + e1-e2;
		sum += math.log(1 + math.pow(10,pvlog-sum),10);
	#print A,R,e,sum;
	return sum;

#print binomial_pvalue(500,40,0.005); sys.exit();

def pvalue_contable_iter(a,stranded,iter,NCRtable):
	cl = compute_likelihood(a);
	pvalue = compute_pvalue_fisher(a,cl,stranded,iter,NCRtable);
	r01 = 0; r00 = 0; 
	for i in range(len(a)): r01 += a[i][0]; r00 += a[i][1];
	fisherpval = 1;
	return [r00,r01,cl,pvalue,fisherpval];

def known_SNPs(snpfile):
	snpindex = {};	
	if not os.path.isfile(snpfile): return snpindex;
	File = open(snpfile,'r');
	for line in File: 	
		snp =line.strip().split();	
		try: 
			vars = int(snp[0]);
			try: 	snpindex[(snp[1],int(snp[2]))] += vars;
			except KeyError: snpindex[(snp[1],int(snp[2]))] = vars;
		except ValueError:
			try:    snpindex[(snp[1],int(snp[2]))] += ':' + snp[0];
                        except KeyError: snpindex[(snp[1],int(snp[2]))] = snp[0];
	File.close(); return snpindex; 

################################################# CRISP MAIN CODE ############################################

bti = {'A':0,'C':1,'G':2,'T':3,'a':5,'c':6,'g':7,'t':8,'N':0,'n':5}; 

QVoffset =33; POOLSIZE = 50; MAXITER = 10000;
MFLAG =2;
readlength = 36;  # works for paired end reads as well.....
MINPOS = 1; MAXPOS =1; MIN_M = 20; MAX_MM = 4;
thresh1 = -3; thresh2 = -5; thresh3 = -2;  # three pvalue thresholds 
MIN_READS = 4; SER = -1;
LOWCOV = 0; 
Qmin=23;# error rate of 0.005

###################################################################################################################


def call_genotypes(ctable,ctable_stranded,ctpval,qvpval,listofsamples,poscounts,ref,alt,poolsize,itb):
	het = float(1)/poolsize; variants =0;
	for a in range(len(ctable)): 
		uniquess = 0; middle = 0;
		if ctable[a][1] > 2:
			for r in range(MINPOS,readlength-MAXPOS): 
				if poscounts[a][r][alt] > 0: uniquess +=1;
				if poscounts[a][r][alt+5] > 0: uniquess +=1;
			for r in range(readlength/4,readlength-readlength/4): 
				if poscounts[a][r][alt] > 0: middle +=1;
				if poscounts[a][r][alt+5] > 0: middle +=1;
			
			pvallowf = minreads_pvalue(ctable_stranded[a][0],ctable_stranded[a][1],min(float(1)/poolsize,0.99));
			pvallowr = minreads_pvalue(ctable_stranded[a][2],ctable_stranded[a][3],min(float(1)/poolsize,0.99));
			minpvlog = pvallowf + pvallowr;
			#if pvallowf > 0 and pvallowr > 0: minpvlog = math.log(pvallowf,10) + math.log(pvallowr,10);
			#else: minpvlog = -100;
			
			pvallowf1 = minreads_pvalue(ctable_stranded[a][0],ctable_stranded[a][1],float(0.75)/poolsize);
			pvallowr1 = minreads_pvalue(ctable_stranded[a][2],ctable_stranded[a][3],float(0.75)/poolsize);
			minpvlog1 = pvallowf1 + pvallowr1;
#			if pvallowf1 > 0 and pvallowr1 > 0: minpvlog1 = math.log(pvallowf1,10) + math.log(pvallowr1,10);
#			else: minpvlog1 = -100;

			# HACK2 for indels
			if (minpvlog >= thresh3 or (minpvlog1 >= thresh3 and ctpval-minpvlog <= thresh1 and poolsize >= 10) or (alt==4 and minpvlog1 >= thresh3 )) and ( (ctable_stranded[a][1] > 0 and ctable_stranded[a][3]> 0) or uniquess >= 6 or alt ==4 or ref ==4)  and ((middle > 0 and uniquess >= MIN_READS) or (alt ==4 and ctable_stranded[a][1] >=1 and ctable_stranded[a][3] >= 1 and ctable[a][1] >= MIN_READS) or ref ==4): variants +=1;

			print '%12s (%d,%d) (%d,%d) %1.2f %1.2f' % (listofsamples[a],ctable_stranded[a][0],ctable_stranded[a][1],ctable_stranded[a][2],ctable_stranded[a][3],pvallowf,pvallowr),
			print '%2d %2d' %(uniquess,middle),
			#print '%2.1f %2.1f %2d %2d' %(qvpval[a][0],qvpval[a][1],uniquess,middle),

	if variants ==1 and ctpval >= -2: variants =0;  # HACK1 if variants = 1 and ctpval is high, ignore 
	print 'variants',variants,
	if ref == 4 or alt ==4: print 'indel',
	else: print 'snp',
	return variants;

def chernoff_bound(counts,stats,refbase,altbase,SER):
	if SER == -1: 
		mu0 = stats[altbase]+stats[refbase]; mu1 = stats[altbase+5]+stats[refbase+5];
		if counts[altbase] ==0 or mu0 ==0: delta0 = 1; 
		else: delta0 = counts[altbase]/mu0-1;
		if counts[altbase+5] ==0 or mu1 ==0: delta1 = 1;
		else: delta1 = counts[altbase+5]/mu1-1;
		if delta0 <= 0: p0 = 0;
		else:	p0 = mu0*(delta0*math.log(math.e,10)-(1+delta0)*math.log(1+delta0,10));
		if delta1 <= 0: p1 = 0;
		else:	p1 = mu1*(delta1*math.log(math.e,10)-(1+delta1)*math.log(1+delta1,10));
		return [p0,p1]; 
	else:
		p0 = binomial_pvalue(counts[refbase]+counts[altbase],counts[altbase],SER);
		p1 = binomial_pvalue(counts[refbase+5]+counts[altbase+5],counts[altbase+5],SER);
		return [p0,p1];
		#if p0 > 0: p0log = math.log(p0,10); 
		#else: p0log = -101;
		#if p1 > 0: p1log = math.log(p1,10); 
		#else: p1log = -101;
		#return [p0log,p1log];
		
		
def indelanalysis(indellist,counts,indcounts): # pairs of the form (sampleid,indelstring) 
	indelcounts = {};
	for i in xrange(len(indellist)):
		snp = indellist[i][1];
		k=1; prev =0;
		while k < len(snp):
			if snp[k] == '+' or snp[k] == '-': 
				if snp[k-1].islower(): 
					try: indelcounts[snp[prev:k].upper()][1] +=1; 
					except KeyError: indelcounts[snp[prev:k].upper()] =[0,1];
				else:
					try: indelcounts[snp[prev:k]][0] +=1; 
					except KeyError: indelcounts[snp[prev:k].upper()] =[1,0];
				prev = k; 
			k +=1;
	
		if snp[k-1].islower(): 
			try: indelcounts[snp[prev:k].upper()][1] +=1; 
			except KeyError: indelcounts[snp[prev:k].upper()] =[0,1];
		else:
			try: indelcounts[snp[prev:k]][0] +=1; 
			except KeyError: indelcounts[snp[prev:k].upper()] =[1,0];

	freqindel = 0; bases= '';
	for indel in indelcounts.iterkeys(): 
		icounts = indelcounts[indel];
		if icounts[0] + icounts[1] > freqindel: freqindel = icounts[0] + icounts[1]; bases = indel;
		#print indel,indelcounts[indel];
	#print 'best',bases,freqindel;
	if freqindel < 3: return ''; # empty string to signify no indel 
	bases0 = bases.lower();
	
	for i in xrange(len(indellist)):
		snp = indellist[i][1];
		k=1; prev =0;
		while k < len(snp):
			if snp[k] == '+' or snp[k] == '-':  
				if snp[prev:k] == bases: indcounts[indellist[i][0]][4] += 1; counts[4] +=1;
				elif snp[prev:k] == bases0: indcounts[indellist[i][0]][9] += 1; counts[9] +=1;
				prev = k;
			k +=1;
		if snp[prev:k] == bases: indcounts[indellist[i][0]][4] += 1; counts[4] +=1; #stats[indellist[i][0]][4] += IER;
		elif snp[prev:k] == bases0: indcounts[indellist[i][0]][9] += 1; counts[9] +=1; #stats[indellist[i][0]][9] += IER;
		#print 'indel',bases,indcounts[indellist[i][0]][4],indcounts[indellist[i][0]][9],indellist[i];

	return bases;

def compute_counts(snp,minQ,MIN_M,MAX_MM,counts,Qhighcounts,indcounts,stats,poscounts):

	sum =0; refbase = snp[2]; rbl = bti[snp[2].lower()]; rbu = bti[snp[2].upper()];
	pvector = snp[7].split(',');

	for t in xrange(1,len(snp[5])):
		if MFLAG ==1: mm = ord(snp[9][t])-33;
		elif MFLAG ==2: mm = ord(snp[8][t-1])-48;   # number of mismatches in read 
		else: mm = 0;
		pos = (int(pvector[t-1])-1)%readlength+1;
		qv = ord(snp[5][t])-QVoffset;

		if qv >= minQ and ord(snp[6][t])-33 >= MIN_M and mm <= MAX_MM  and pos > MINPOS and pos <= readlength-MAXPOS: 
			if snp[4][t] == '.' and mm <= MAX_MM-1: 
				counts[rbl] +=1; indcounts[rbl] += 1;stats[rbl] += math.pow(0.1,0.1*qv);
				if qv >= Qmin: Qhighcounts[rbl] +=1;
				poscounts[pos][rbl] +=1;
			elif snp[4][t] == ',' and mm <= MAX_MM-1: 
				counts[rbu] +=1; indcounts[rbu] += 1;stats[rbu] += math.pow(0.1,0.1*qv);
				if qv >= Qmin: Qhighcounts[rbu] +=1;
				poscounts[pos][rbu] +=1;
			elif snp[4][t] != '.' and snp[4][t] != ',': 
				try:
					base = bti[snp[4][t]]; 
					if qv >= Qmin: Qhighcounts[base] +=1;
					counts[base] +=1; indcounts[base] +=1;  stats[base] += math.pow(0.1,0.1*qv);
					poscounts[pos][base] +=1;
				except KeyError: pass;
			sum +=1;
	#if len(snp) >= 10: pass; # concatenate to a string of indels separated by ':' 

#	print snp[0],snp[1],snp[2],sum,indcounts; 
	return sum;

def compare_counts_pools(dirname,snpfile,poolsize,minQ):
	maxN = 1000;  NCRtable = compute_NCRtable(maxN);
	print >>sys.stderr, 'finished computing NCR table';
	MINQpv = thresh2;

	snpindex = known_SNPs(snpfile);
	files = []; filehandlers = [];	listofsamples = []; 

	csl = dirname.split(',');
	if len(csl) > 1: 
		for file in csl:
			files.append(os.path.basename(file)); filehandlers.append(open(file,'r'));  
			samplename = os.path.basename(file).rstrip('.pileup'); listofsamples.append(samplename);
		noofsamples = len(csl);
			
	else:

		if os.path.isdir(dirname): dirname += '/*pileup';
		for file in glob.glob(dirname): 
			files.append(os.path.basename(file)); filehandlers.append(open(file,'r')); 
			samplename = os.path.basename(file).rstrip('.pileup');
			#samplename = os.path.basename(file).split('.')[0] + '.' + os.path.basename(file).split('.')[1];
			listofsamples.append(samplename);
		noofsamples = len(files);


	print >>sys.stderr, "noofsamples:",noofsamples; #sys.exit();
	if noofsamples < 2: print >>sys.stderr, 'CRISP requires at least 2 pools of sequence data (pileup files) to make variant calls'; sys.exit();

	#print listofsamples;
	itb = {0:'A',1:'C',2:'G',3:'T'};
	finish = 0;
	poscounts = []; indcounts = [];  counts = [0,0,0,0,0,0,0,0,0,0]; stats = []; tstats = []; Qhighcounts  = [0,0,0,0,0,0,0,0,0,0]; 
	for q in range(10): tstats.append(0.0);
	for i in range(noofsamples):
		stats.append([]);
		poscounts.append([]); indcounts.append([0,0,0,0,0,0,0,0,0,0]);
		for q in range(readlength+1): poscounts[i].append([0,0,0,0,0,0,0,0,0,0]);
		for q in range(10): stats[i].append(0.0);

	print '#non 0 CHROM\tposition\trefbase\tallele1\tallele2\t allele2-count+ coverage+ allele2-count- coverage- allele3 allele3-count | CT CT-pvalue CHI2 chisq-pvalue QV QVpvalue+ QVpvalue- SD strand-bias';
	positions = 0;
	while finish ==0:
		
		refb = 'N'; qpvalue = []; qpvars =0; qpsum =[0.0,0.0];
		for q in range(10): counts[q] = 0; tstats[q] = 0.0; Qhighcounts[q] =0;
		indellist = []; iflag =0;

		for i in range(noofsamples): 
			newline = filehandlers[i].readline(); 
			if not newline: finish = 1;  break;
			snp = newline.split(); 
			for q in range(10):  indcounts[i][q] = 0; stats[i][q] = 0.0;
			for q in range(readlength+1):
				for p in range(10):  poscounts[i][q][p] = 0;
			if snp[2] != 'N': refb = snp[2];
			if int(snp[3]) > 0: sum = compute_counts(snp,minQ,MIN_M,MAX_MM,counts,Qhighcounts,indcounts[i],stats[i],poscounts[i]);
			else: sum =0;
			if len(snp) >= 10: indellist.append([i,snp[9]]); iflag += len(snp[9])/2;
		if finish == 1: break;
		else: positions +=1;
		if positions%100000 ==0: print >>sys.stderr, 'processed',positions,'bases';

		idallele = '';
		if iflag  >= 3: 
			idallele = indelanalysis(indellist,counts,indcounts);
			Qhighcounts[4] = counts[4]; Qhighcounts[9] = counts[9]; # for indel, copy counts


		maxvec = [];
#		for base in ['A','C','G','T']:	maxvec.append([counts[bti[base]]+counts[bti[base]+4],base]);
		for i in [0,1,2,3]: maxvec.append([counts[i]+counts[i+5],i]);
		#for i in [0,1,2,3]: maxvec.append([Qhighcounts[i]+Qhighcounts[i+5],i]);
		if idallele != '':  maxvec.append([counts[4]+counts[4+5],4]); 	itb[4] = idallele;
		#print maxvec;

		maxvec.sort(); maxvec.reverse(); allele1 = maxvec[0][1]; allele2 = maxvec[1][1]; allele3 = maxvec[2][1];
		refbase = bti[refb];

		asymptotic_table =1; 
		ctable = [];		ctable_stranded = []; ctablef = []; ctabler = []; 
		for i in range(noofsamples): 
			c0 = indcounts[i][allele1] + indcounts[i][5+allele1];
			c1 = indcounts[i][allele2] + indcounts[i][5+allele2];
			if c0 < 5 or c1 < 5: asymptotic_table = 0;
			if allele1 == refbase: # c1 is always alternate base  
				ctable.append([c0+c1,c1]);
				ctable_stranded.append([indcounts[i][allele1] + indcounts[i][allele2],indcounts[i][allele2],indcounts[i][5+allele1] + indcounts[i][5+allele2],indcounts[i][5+allele2]]);
			else: 
				ctable.append([c0+c1,c0]);
				ctable_stranded.append([indcounts[i][allele1] + indcounts[i][allele2],indcounts[i][allele1],indcounts[i][5+allele1] + indcounts[i][5+allele2],indcounts[i][5+allele1]]);
		if (snp[0],int(snp[1])) in snpindex: print 'snp',snpindex[(snp[0],int(snp[1]))],
		else: print 'non %2d' % (0),
		print '%6s %6s %1s %1s %1s ' % (snp[0],snp[1],snp[2],itb[allele1],itb[allele2]), 

		print '%5d %5d %5d %5d' % (counts[allele2],counts[allele1]+counts[allele2],counts[allele2+5],counts[allele1+5]+counts[allele2+5]),
		print '%1s %3d |' %(itb[allele3],counts[allele3]+counts[allele3+5]),

		if allele1 != refbase: 
			altbase = allele1;
			if allele2 != refbase: refbase = allele2;  # special case, refbase is not seen 
		else: altbase =  allele2;

		if altbase ==4 or refbase ==4: ser = 0.01; # HAVE to use quality scores for indels
		else: ser = SER; 
		if LOWCOV ==0: 
			for i in range(noofsamples): 
				qpvalue.append(chernoff_bound(indcounts[i],stats[i],refbase,altbase,ser));
				if qpvalue[i][0] <= MINQpv and qpvalue[i][1] <= MINQpv: qpvars +=1; qpsum[0] += qpvalue[i][0]; qpsum[1] += qpvalue[i][1];
			pvalue_all = chernoff_bound(Qhighcounts,stats,refbase,altbase,math.pow(0.1,0.1*Qmin));

		#if maxvec[1][0] >= 3 or (allele1 != refbase and maxvec[0][0] >= 3):  changed june 16 2010
		if counts[allele2] + counts[allele2 +5] >= 3 or (allele1 != refbase and counts[allele1] + counts[allele1 +5] >=3): 
			r00 = 100; r01 = 0; cl = -1; pvalue = 0.1; pvalue_allreads = 0.1; strandpval = 0.1;
			if len(ctable) ==2: 
				#[r00,r01,cl,pvalue,pvalue_allreads] = pvalue_contable_iter(ctable,0,10,NCRtable);
				pvalue = fet(ctable[0][1],ctable[0][0],ctable[1][1],ctable[1][0]);
				#print pvalue,ctable;
			else: 
				if asymptotic_table ==1: [r00,r01,cl,pvalue,pvalue_allreads] = pvalue_contable_iter(ctable,0,1000,NCRtable);
				else: [r00,r01,cl,pvalue,pvalue_allreads] = pvalue_contable_iter(ctable,0,MAXITER,NCRtable);

			strandpval = fet(counts[allele2],counts[allele2]+counts[allele1],counts[allele2+5],counts[allele2+5]+counts[allele1+5]);
			pvaluelog = -50; 
			if pvalue > 0: pvaluelog = math.log(pvalue,10);
			pvalchi2 = pvalue_chi2(ctable);
			if pvalchi2 > 0: pvalchi2log = math.log(pvalchi2,10);
			else: pvalchi2log = -50;
			if asymptotic_table ==1: jointpvallog = pvalchi2log;
			else: jointpvallog = pvaluelog;
			if jointpvallog  <= thresh1 or (qpvars >= 1) or (pvalue_all[0] <= -10 and pvalue_all[1] <= -10 and pvalue <= 0.05):
				#if asymptotic_table ==1 and -1*pvaluelog <= math.log(MAXITER,10) : jointpvallog = pvalchi2log;
				strandpvallog = 0;
				if strandpval > 0:  strandpvallog = math.log(strandpval,10);
				print '  CT %+2.1f chisq-%d %+2.1f QV %+2.1f %+2.1f %2.1f %2.1f SD %+2.1f |' % (pvaluelog,asymptotic_table,pvalchi2log,qpsum[0],qpsum[1],pvalue_all[0],pvalue_all[1],strandpvallog),
				#for i in xrange(len(ctable)): print '%2.1f %2.1f ' %(qpvalue[i][0],qpvalue[i][1]),
				variants = call_genotypes(ctable,ctable_stranded,jointpvallog,qpvalue,listofsamples,poscounts,refbase,altbase,poolsize,itb);
				if variants > 0: 
					print '| allelecounts', 
					for i in range(noofsamples): 
						c0 = indcounts[i][allele1] + indcounts[i][5+allele1];
						c1 = indcounts[i][allele2] + indcounts[i][5+allele2];
						print listofsamples[i] + ':' + `c0` + ':' + `c1`,
				print;
			else: print '  CT %+2.1f QV %+2.1f %+2.1f %2.1f %2.1f ' %(math.log(pvalue,10),qpsum[0],qpsum[1],pvalue_all[0],pvalue_all[1]);
		else: print;

	for i in range(noofsamples): filehandlers[i].close();


#################################################################################################################################

parser = OptionParser();
parser.add_option("-d","--pileupdir",dest="pileupdir",help="directory with pileup files");
parser.add_option("--rl","--readlength",type="int",dest="readlength",help="length of reads",default=36);
parser.add_option("--qvoffset",dest="Qoffset",type="int",help="quality value offset, 33 for Sanger, 64 for Illumina base qualities",default=33);
parser.add_option("--mbq",dest="MIN_Q",type="int",help="minimum base quality to consider a base for snp calling, default 10",default=10);
parser.add_option("--mmq",dest="MIN_M",type="int",help="minimum read mapping quality to consider a read for snp calling, default 20",default=20);
parser.add_option("--maxm",dest="MAX_MM",type="int",help="maximum number of mismatches allowed for read to be considered for snp calling. The number of mismatches for each read is stored in column 9 of the pileup file. Default value is max(3,0.04 x readlength)",default=3);
parser.add_option("--dbfile",dest="ncbifile",help="list of dbSNP variants in sequenced region",default="missing");
parser.add_option("--clipl",dest="MIN_POS",type="int",help="ignore the first x bases of a read for snp calling, default 0",default=0);
parser.add_option("--clipr",dest="MAX_POS",type="int",help="ignore the last x bases of a read for snp calling, default 0",default=0);
parser.add_option("--minc",dest="MIN_READS",type="int",help="minimum number of reads with alternate allele required for calling position as SNP, default 4",default=4);
parser.add_option("-p","--poolsize",dest="POOLSIZE",type="int",help="number of haplotypes in each pool",default=2);
parser.add_option("--ctpval",dest="thresh1",type="float",help="threshold on the contingency table p-value for calling position as variant (specified as log10)",default=-3);
parser.add_option("--qvpval",dest="thresh2",type="float",help="threshold on the quality values p-value for calling position as variant (specified as log10)",default=-5);
parser.add_option("--perms",dest="MAXITER",type="int",help="maximum number of permutations for contingency table p-value, default 10K",default=10000);
parser.add_option("-e","--seqerate",dest="SER",type="float",help="estimate of average sequencing error rate, used for quality value based p-value. Default value is -1 which uses individual base qualities to estimate p-values",default=-1);

(options,args) = parser.parse_args();
MIN_Q = options.MIN_Q; MIN_M = options.MIN_M; 
MINPOS = options.MIN_POS; MAXPOS = options.MAX_POS; 
MIN_READS = options.MIN_READS;  POOLSIZE = options.POOLSIZE; MAXITER = options.MAXITER;
readlength = options.readlength; QVoffset = options.Qoffset;
MAX_MM = max(options.MAX_MM,int(0.04*readlength));
SER = options.SER;
thresh1 = options.thresh1;thresh2 = options.thresh2;

print >>sys.stderr, 'options: readlength',readlength,' haplotypes_per_pool',POOLSIZE,'min_base_quality',MIN_Q,'max_mismatches',MAX_MM,'CT_pvalue_max',thresh1,'qvalue_pvalue_max',thresh2,'quality-value-offset',QVoffset;

if options.pileupdir == None:
	print '\n-----------------------------------------------------------------------------------------------------';
	print '                      python CRISP.py -h ';
	print '\n-----------------------------------------------------------------------------------------------------'
	sys.exit();

else: compare_counts_pools(options.pileupdir,options.ncbifile,POOLSIZE,MIN_Q);

