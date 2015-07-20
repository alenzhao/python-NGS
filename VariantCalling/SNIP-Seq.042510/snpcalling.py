#!/usr/bin/python2.4
import sys, os, glob, string, subprocess,time, math, re, compiler, random
space = re.compile(r'\s+');

pflag = 0;
MFLAG = 1;
#Qoffset=33;
Moffset=33;
MMoffset=48;
LOW_MQ = 0; 

btimap = {',':[-1,0], '.':[-1,4],'A':[0,0],'C':[1,0],'G':[2,0],'T':[3,0],'a':[0,4],'c':[1,4],'g':[2,4],'t':[3,4],'n':[-2,0],'N':[-2,0]};
qsmap = [];
for t in range(10): qsmap.append(0);
for t in range(10): qsmap.append(1);
for t in range(10): qsmap.append(2);
for t in range(40): qsmap.append(3);

def bti(b): # base to int 
	if b == 'A' or b == 'a': return 0;
	elif b == 'C' or b == 'c': return 1;
	elif b == 'G' or b == 'g': return 2;
	elif b == 'T' or b == 't': return 3;
	else: return -1;

def itb(b): # int to base
	if b == 0: return 'A';
	elif b == 1: return 'C';
	elif b == 2: return 'G';
	elif b == 3: return 'T';
	else: return 'N';


def filter_bin(genposition,refbase,altbase,position,strand,UniqueTable,RL): # decide if we should filter this bin based on UniqueTable
	return 0;
	if strand == 4: position = RL - position-1; 
	if len(UniqueTable[genposition]) == 0: return 0; 
	if UniqueTable[genposition][position][altbase] ==1: return 1; 
	if UniqueTable[genposition][position][altbase] ==2: return 2; 
	else: return 0; 

def print_bintable_vertical(qcounts,posbins):
	for r in range(posbins):
		print 'F %2d |' % (r+1),
		for k in range(4):
			for t in range(4): print '%3d' % (qcounts[r*8+t][k]),
			print '|',
		print '-- R %2d |' % (r+1),
		for k in range(4):
			for t in range(4): print '%3d' % (qcounts[r*8+t+4][k]),
			print '|',
		print;


def compute_qvalues(snp,qcounts,counts,qsums,j,posbins,potsnplist,indcounts,refbase,MAQSNP,Qoffset):
	pvector = snp[7].split(','); recal_qscores = '@';
	for k in range(1,len(snp[4])):
		ps = (int(pvector[k-1])-1); 
		quality = ord(snp[5][k])-Qoffset; 	qoff = qsmap[quality]; 
#		if ord(snp[6][k])-Moffset < LOW_MQ: qoff = 0; # added to handle low mapping quality reads for this bin  
		bquality = pow(0.1,float(quality)/10); 
		[base,offset] = btimap[snp[4][k]]; 
		if base == -1: base = refbase; 
		empirical_erate = 0; deltasum = 0; deltaq = 0;

		# initial prior = 1,x where x is quality value NEW CODE FOR beta prior added juy 16 2009 
		initialprior = int(float(1)/bquality);
		initbeta = 20;
		initalpha = float(bquality)*float(initbeta);
		if initalpha ==0: initalpha = float(initbeta)/100;
		maxCOV = 0.01; delta = 0.01;
		if base != refbase:
			betamean = float(qcounts[8*ps+base+offset][qoff]+ delta)/(qsums[2*ps+offset/4][qoff]+delta*2); 
			weight = math.sqrt(qsums[2*ps+offset/4][qoff])/7;  #max weight set to 50
			if weight > 1: weight = 1;
			newerrorrate = bquality + (1-bquality)*weight*betamean; 
			new_qscore = math.log(newerrorrate,10)*(-10); recal_qscores += chr(int(new_qscore)+Qoffset);
			
		else:
			betamean = float(qsums[2*ps+offset/4][qoff]-qcounts[8*ps+base+offset][qoff]+ delta)/(qsums[2*ps+offset/4][qoff]+delta*2); 
			weight = math.sqrt(qsums[2*ps+offset/4][qoff])/7;
			if weight > 1: weight = 1;
			newerrorrate = bquality + (1-bquality)*weight*betamean; 
			new_qscore = math.log(newerrorrate,10)*(-10); recal_qscores += chr(int(new_qscore)+Qoffset);

	return recal_qscores; 

############################################################################################################################3

# function for computing counts of A/C/T/G in each bin 

def filter_snps(snp):
	alts =0;
	for k in range(1,len(snp[4])):
		[base,offset] = btimap[snp[4][k]]; 
		if base != -1: alts +=1;
	if alts*10 >= len(snp[4])-1: return alts;
	else: return 0;

def compute_counts(snp,counts,qsums,j,refbase,indcounts,posbins,filename,qcounts,var,MIN_Q,MIN_M,MAX_MM,RL,Qoffset):
	pvector = snp[7].split(','); 
	for k in range(1,len(snp[4])):
		ps = (int(pvector[k-1])-1); 
#		if ord(snp[5][k]) > 73 and Qoffset == 33: Qoffset = 64; 
		qoff = qsmap[ord(snp[5][k])-Qoffset]; 
#		if ord(snp[6][k])-Moffset < LOW_MQ: qoff = 0; # added to handle low mapping quality reads for this bin  
		[base,offset] = btimap[snp[4][k]]; 
		if base == -1: base = refbase; 
		if ord(snp[6][k])- Moffset < MIN_M: continue; # mapping score 
		if ord(snp[5][k])- Qoffset < MIN_Q: continue; # quality score
		if MAX_MM > 0: mm = ord(snp[8][k-1])-MMoffset; # if MAX_MM ==0, column 9 is missing  
		else: mm = -2;
		if RL  == 0: MAX_MM = ord(snp[9][k-1])/17 +1;
		if base == refbase and mm > MAX_MM-1: continue;
		if base != refbase and mm > MAX_MM: continue;
		counts[8*(posbins)+base+offset][qoff] +=1;  counts[8*ps+ base+offset][qoff] += 1; 
		indcounts[j][8*(posbins)+base+offset][qoff] +=1;  indcounts[j][8*ps+base+offset][qoff] +=1;

	# determine two most likely bases based on Q20 and Q10 counts alone
	cvector = [];
	for k in range(4): cvector.append([indcounts[j][8*posbins+k][2]+indcounts[j][8*posbins+k+4][2]+indcounts[j][8*posbins+k][1]+indcounts[j][8*posbins+k+4][1],k]); 
	cvector.sort(); 
	if cvector[3][1] != refbase:
		if cvector[3][0] >= 3: return [1,cvector[3][1],len(pvector)];
	elif (cvector[2][0] >= 2 and cvector[2][0] >= 0.2*cvector[3][0]) or cvector[2][0] >= 3: return [1,cvector[2][1],len(pvector)];
	return [0,-1,len(pvector)]; 


# compute three genotype probabilities for each individual and a position 'j' for each strand separately 
def compute_genotypes(snp,samplename,printflag,ref,alt,het,delta,relax,qvector,MIN_Q,MIN_M,MIN_POS,MAX_POS,MAX_MM,posbins,readstosample,UniqueTable,genpos,RL,Qoffset):
	bvector = snp[4]; mvector = snp[6]; pvector = snp[7].split(',');
	refbase = bti(snp[2]); position = int(snp[1]);
	counts = [0,0,0,0,0,0,0,0]; mscores = [0,0,0,0,0,0,0,0]; 	qscores = [0,0,0,0,0,0,0,0];
	mmcounts = [0,0,0,0,0,0,0,0];	lowmapcounts = [0,0,0,0,0,0,0,0]; 	read2counts = [0,0,0,0,0,0,0,0];	pcounts = [0,0,0,0,0,0,0,0]; # position counts 
	posdist = [];
	for i in range(8): posdist.append([0,0,0]); 
	sample = os.path.basename(samplename).split('.')[0] + '.'+ os.path.basename(samplename).split('.')[1];
	maxMQ = [0,0,0,0,0,0,0,0]; 

	reads =0; # no of filtered reads 
	for k in range(1,len(bvector)):
		[base,offset] = btimap[snp[4][k]];
		if base == -1: base = refbase; 
		if ord(snp[5][k])-Qoffset < MIN_Q: continue;  # this represents the original quality 
		if ord(snp[6][k])-Moffset < MIN_M: lowmapcounts[base+offset] +=1; continue;
		if MAX_MM > 0: mm = ord(snp[8][k-1])-MMoffset; 
		else: mm = -2;
		if RL == 0: MAX_MM = ord(snp[9][k-1])/17 +1; readlength = ord(snp[9][k-1]);	
		else: readlength = RL;
		p = (int(pvector[k-1])-1)%readlength +1;
		if p <= MIN_POS or readlength -p  < MAX_POS: continue;  # position filter
		if base == refbase and mm > MAX_MM-1: continue;
		if base != refbase and mm > MAX_MM: continue;
		if ord(snp[6][k])-Moffset > maxMQ[base+offset]: maxMQ[base+offset] = ord(snp[6][k])-Moffset; 
		reads +=1; 
		counts[base+offset] += 1; mscores[base+offset] += ord(mvector[k])-Moffset; qscores[base+offset] += ord(snp[5][k])-Qoffset; 
		if int(pvector[k-1])%readlength < 8: posdist[base+offset][0] +=1;
		elif int(pvector[k-1])%readlength >= readlength-8: posdist[base+offset][2] +=1;
		else: posdist[base+offset][1] +=1; 
		pcounts[base+offset] += min(readlength-int(pvector[k-1])%readlength,int(pvector[k-1])%readlength); # distance from end of read 
		#if snp[8][k-1] == '2': read2counts[base+offset] +=1; 
# average position, average # of mismatches, average counts on read1 or read2 

	[ref,alt] = besttwo_qscores(qscores,refbase); # use 0 instead of posbins for making it work
	plist = {}; uniquepositions0 = 0; uniquepositions1 = 0;
	g00f = 0; g01f = 0; g11f =0; g00r =0; g01r =0; g11r=0;
	ERRf=0.0; ERRr=0; ERAf=0; ERAr=0;
	for k in range(1,len(bvector)):
		[base,offset] = btimap[snp[4][k]];
		if base == -1: base = refbase; 
		if ord(snp[5][k])- Qoffset < MIN_Q or ord(snp[6][k])- Moffset < MIN_M: continue; # mapping score 
		if MAX_MM > 0: mm = ord(snp[8][k-1])-MMoffset; 
		else: mm = -2;
		if RL == 0: MAX_MM = ord(snp[9][k-1])/17 +1; readlength = ord(snp[9][k-1]);
		else: readlength = RL;
		p = (int(pvector[k-1])-1)%readlength +1;
		if p <= MIN_POS or readlength -p  < MAX_POS: continue;  # position filter
		if base == refbase and mm > MAX_MM-1: continue;
		if base != refbase and mm > MAX_MM: continue;

		if bvector[k] == ',': 
			p = pow(10,float(Qoffset-ord(qvector[k]))/10); g01f += math.log(1-het,10);
			if p < 1: g00f += math.log(1-p,10); 
	#		g11f += float(Qoffset-ord(snp[5][k]))/10;
			g11f += float(Qoffset-ord(qvector[k]))/10;  # should this be changed to snp[k] since we dont estimate the probability of reading alternate base as reference base 
		elif bvector[k] == '.': 
			p = pow(10,float(Qoffset-ord(qvector[k]))/10); g01r += math.log(1-het,10);
			if p < 1: g00r += math.log(1-p,10); 
			g11r += float(Qoffset-ord(qvector[k]))/10; 
		elif offset ==0 and base ==alt: 
			p = pow(10,float(Qoffset-ord(qvector[k]))/10); g00f += float(Qoffset-ord(qvector[k]))/10; g01f += math.log(het,10);
			#if snp[8][k-1] == '2': pos = int(pvector[k-1])+36; 
			#else: pos = int(pvector[k-1]); 
			pos = int(pvector[k-1]); 
			if (pos,0) not in plist: uniquepositions0 +=1; plist[(pos,0)] = 1; 
			if p < 1: g11f += math.log(1-p,10); 
		elif offset ==4 and base ==alt: 
			p = pow(10,float(Qoffset-ord(qvector[k]))/10); g00r += float(Qoffset-ord(qvector[k]))/10; g01r += math.log(het,10);
			if p < 1: g11r += math.log(1-p,10);
			#if snp[8][k-1] == '2': pos = int(pvector[k-1])+36; 
			#else: pos = int(pvector[k-1]); 
			pos = int(pvector[k-1]); 
			if (pos,1) not in plist: uniquepositions1 +=1; plist[(pos,1)] = 1; 

	max = 0;
	if g01f + g01r >= g00f + g00r + delta and g01f >= g00f + relax and g01r >= g00r + relax: max = 1;
	# additional loose filter for SNPs 
	delta = 6; relax = -1;  
	if g01f + g01r >= g00f + g00r + delta and g01f >= g00f + relax and g01r >= g00r + relax: max = 1;
	if max ==1 and g11f + g11r >= g01f + g01r + 1: max = 2; # call heterozygote as homozygote
	
	#if maxMQ[alt] < 30 or maxMQ[alt+4] < 30: max = 0; # filter for low max mapping quality added July 5 2009 
	
	#if sample in snp_known and printflag ==1:  printflag = 2; print 'maq 1',
	if (max >=1 and printflag ==1) or (printflag ==2): 
		#if sample not in snp_known: print 'maq 0',
		print itb(ref),itb(alt),itb(refbase), 
		print 'ICALL %4d %20s %1s %1s' % (position,os.path.basename(sample).rstrip('.pileup'),itb(ref),itb(alt)),
		for t in range(4): print '%2d %2d' % (counts[t],counts[t+4]),
		print 'UNIQUE %2d %2d' % (uniquepositions0,uniquepositions1), 
		print 'MS',
		for t in [ref,alt]: print '%2.1f %2.1f' % (float(mscores[t])/(counts[t]+0.01),float(mscores[t+4])/(counts[t+4]+0.01)),
		print 'MAXMQ',
		for t in [ref,alt]: print '%2d %2d' % (maxMQ[t],maxMQ[t+4]),
		print 'BQ',
		for t in [ref,alt]: print '%2.1f %2.1f' % (float(qscores[t])/(counts[t]+0.01),float(qscores[t+4])/(counts[t+4]+0.01)),
		print 'indel',
		for t in [ref,alt]: print '%2.1f %2.1f' % (float(pcounts[t])/(counts[t]+0.01),float(pcounts[t+4])/(counts[t+4]+0.01)),
		print 'MISMATCH',
		print '%2.1f' % (float(mmcounts[alt])/(counts[alt]+0.01)-1-float(mmcounts[ref])/(counts[ref]+0.01)),
		print '%2.1f' % (float(mmcounts[alt+4])/(counts[alt+4]+0.01)-1-float(mmcounts[ref+4])/(counts[ref+4]+0.01)),
		#print 'Lowmap',
		#for t in [ref,alt]: print '%2.1f %2.1f' % (float(lowmapcounts[t])/(counts[t]+lowmapcounts[t]+0.01),float(lowmapcounts[t+4])/(counts[t+4]+lowmapcounts[t+4]+0.01)),
		print 'read2',
		for t in [ref,alt]: print '%2.1f %2.1f' % (float(read2counts[t])/(counts[t]+0.01),float(read2counts[t+4])/(counts[t+4]+0.01)),
		print '  PROB','%2.1f %2.1f %2.1f %2.1f %2.1f %2.1f ' % (g00f,g01f,g11f,g00r,g01r,g11r),
		print 'PHRED %2d %2d ' % (int(g01f*10-g00f*10),int(g01r*10-g00r*10)),
		if max ==1: print 'HET'; 
		elif max == 2: print 'ALT';

		for k in range(1,len(bvector)):
			[base,offset] = btimap[snp[4][k]];
			if base == -1: base = refbase; 
			ps = int(pvector[k-1])-1; 
			if ord(snp[5][k])-Qoffset < MIN_Q or ord(snp[6][k])-Moffset < MIN_M: continue; 
			if MAX_MM > 0: mm = ord(snp[8][k-1])-MMoffset; 
			else: mm = -2;
			if RL ==0: MAX_MM = ord(snp[9][k-1])/17 +1; readlength = ord(snp[9][k-1]);
			else: readlength = RL;
			p = (int(pvector[k-1])-1)%readlength +1;
			if p <= MIN_POS or readlength -p  < MAX_POS: continue;  # position filter
			if base == refbase and mm > MAX_MM-1: continue;
			if base != refbase and mm > MAX_MM: continue;
			if bti(snp[4][k]) == alt and counts[alt]+counts[alt+4]<15: 
				if filter_bin(genpos,refbase,alt,int(pvector[k-1])-1,offset,UniqueTable,RL) > 0: print 'FB',
				print snp[4][k],ord(snp[5][k])-Qoffset,ord(qvector[k])-Qoffset,'MS',ord(snp[6][k])-Moffset,ps,'|',
	#sys.stderr.read();
		print;
	return [max,ref,alt,g01f-g00f,g01r-g00r,g11f-g01f,g11r-g01r,counts,[uniquepositions0,uniquepositions1,maxMQ[alt],maxMQ[alt+4],int(pcounts[alt]/(counts[alt]+0.01)),int(pcounts[alt+4]/(counts[alt+4]+0.01)),posdist[alt],posdist[alt+4]]];

###################################################################################################################

def Cnr(n,r):
	comb =0.0;
	if 2*r > n: return Cnr(n,n-r); 
	for i in range(r):
		comb += math.log(n-i,10);
		comb -= math.log(r-i,10);
	return comb;

def besttwo_qscores(qscores,refbase):
	scorelist = [];
	for k in range(4): 	scorelist.append([qscores[k] + qscores[4+k],k]);
	scorelist.sort();
	if scorelist[3][1] == refbase: return [scorelist[3][1],scorelist[2][1]]; 
	if scorelist[2][1] == refbase: return [scorelist[2][1],scorelist[3][1]]; 
	if scorelist[1][1] == refbase: return [scorelist[1][1],scorelist[3][1]]; 
	if scorelist[0][1] == refbase: return [scorelist[0][1],scorelist[3][1]]; 

def besttwoQ20bases(counts,posbins):
	scorelist = [];
	for k in range(4): 	scorelist.append([counts[8*posbins+k][2] + counts[8*posbins+4+k][2],k]);
	scorelist.sort();
	return [scorelist[3][1],scorelist[2][1]]; 

def besttwo(counts,refbase,posbins):
	scorelist = [];
	for k in range(4): 	scorelist.append([counts[8*posbins+k] + counts[8*posbins+4+k],k]);
	scorelist.sort();
	return [scorelist[3][1],scorelist[2][1]]; 

###################################################################################################################


def flankingseq(refsequence,i):
	s1 = max(0,i-20); s2 = min(i+20,len(refsequence)-1); seq = '' + refsequence[s1:i] + ' ' + refsequence[i] + ' ' + refsequence[i+1:s2+1];
	return seq;


def read_MAQsnpfile(SNPfile,SNPlist_MAQsnps):
	if not os.path.isfile(SNPfile):  return -1;
	File = open(SNPfile,'r'); MAQsnps =0;
	for f in File:  
		snp = f.split(); 
		sample = os.path.basename(snp[0]).split('.')[0] + '.'+os.path.basename(snp[0]).split('.')[1];
		try: SNPlist_MAQsnps[snp[1]+':'+snp[2]] +=1; 
		except KeyError: SNPlist_MAQsnps[snp[1]+':'+snp[2]] = 1; MAQsnps +=1; 
		SNPlist_MAQsnps[snp[1] + ':'+ snp[2] + ':' + sample] = snp[3:]; 
	File.close(); 
	return MAQsnps;

def read_dbsnpfile(ncbifile,NCBIlist):
	if not os.path.isfile(ncbifile): return -1;
	File = open(ncbifile,'r');
	for f in File:  snp = f.split(); NCBIlist[snp[1]+':'+snp[2]]= snp[4]+'/'+snp[0];
	File.close();
###################################################################################################################



def init_arrays(counts,posbins,indcounts,samples):
	for r in range(8*posbins+8): counts[r][0] = counts[r][1] = counts[r][2] = counts[r][3] = 0;
	for j in range(samples):
		for r in range(8*posbins+8): indcounts[j][r][0] = indcounts[j][r][1] = indcounts[j][r][2] = indcounts[j][r][3] =  0;



def bestMAQsnp(postolocus,filenamelist,SNPlist_MAQsnps):
	maqsnps = 0; bestPHRED = 0;
	for j in range(len(filenamelist)):
		sample = os.path.basename(filenamelist[j]).split('.')[0] + '.'+ os.path.basename(filenamelist[j]).split('.')[1];
		try: 
			snpinfo = SNPlist_MAQsnps[postolocus[0]+ ':'+ postolocus[1]+ ':'+ sample];
			if int(snpinfo[2]) > bestPHRED: bestPHRED = int(snpinfo[2]); 
			maqsnps += 1;
		except KeyError: pass; 
	return [maqsnps,bestPHRED];


def ncr(n,r):
	k = min(r,n-r); ll = 0;
	for i in range(k): ll += math.log(float(n-i)/(i+1),2);
	return ll;

def HWE_exact(refs,hets,alts):
	N = refs + hets + alts;
	n0 = alts*2 + hets;
	n01 = hets; 
# N is sample size, n0 is number of minor allele, n01 is number of heterozygotes 
# Plow = P(N01 <= n01 | N,n0)
	if n0 > N: n0 = 2*N-n0; 
	p = float(n0)/(2*N);	plow = 0; phigh =0; constant = ncr(2*N,n0);
	if n0%2 ==1: # odd
		for N01  in range(1,n01+1,2):
			N0 = n0; N00 = (n0-N01)/2; N11 = N - N01 - N00; N1 = 2*N-N0; 
			plow += math.pow(2,N01 + ncr(N,N00) + ncr(N-N00,N01) - constant); 
		for N01  in range(n01,n0+1,2):
			N0 = n0; N00 = (n0-N01)/2; N11 = N - N01 - N00; N1 = 2*N-N0; 
			phigh += math.pow(2,N01 + ncr(N,N00) + ncr(N-N00,N01) - constant); 

	if n0%2 ==0: # even
		for N01  in range(0,n01+1,2):
			N0 = n0; N00 = (n0-N01)/2; N11 = N - N01 - N00; N1 = 2*N-N0; 
			plow += math.pow(2,N01 + ncr(N,N00) + ncr(N-N00,N01) - constant); 
		for N01  in range(n01,n0+1,2):
			N0 = n0; N00 = (n0-N01)/2; N11 = N - N01 - N00; N1 = 2*N-N0; 
			phigh += math.pow(2,N01 + ncr(N,N00) + ncr(N-N00,N01) - constant); 

#	print 2*p*(1-p)*N,plow,phigh;
	return [2*p*(1-p)*N,plow,phigh];


def compute_HWE(SNPlist):
	countAA = 0; countAa = 0; countaa = 0;
	for i in range(len(SNPlist)):
		if SNPlist[i][0] == 0 : countAA +=1;
		if SNPlist[i][0] == 1 : countAa +=1;
		if SNPlist[i][0] == 2 : countaa +=1;
	[expectedhet,plow,phigh] = HWE_exact(countAA,countAa,countaa);
	return [countAA,countAa,countaa,math.log(min(plow,phigh),10)];


#HWE_exact(165,2,4);
