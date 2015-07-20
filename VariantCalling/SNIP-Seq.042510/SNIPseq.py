#! /usr/bin/env python
# AUTHOR VIKAS BANSAL last edited Feb 3 2010
# SNIP-Seq: Program TO CALL SNPS from POPULATION DATA using multiple PILEUP FILES 

import sys, os, glob, string, subprocess,time, math, re, compiler, random
from optparse import OptionParser
from array import *
import snpcalling
space = re.compile(r'\s+'); compiler.parseFile(sys.argv[0]); 
printflag = 0;
#posbins = 36;   #maximum number of read cycles, i.e. length of read, this is an input argument 
hetratio = 0.5; # fraction of reads with alternate allele for heterozygote..
#MIN_Q = 10; # don't consider base-calls with quality below this value 
#MIN_M = 10; # don't consider reads with mapping score below this value
#MIN_POS = 2; MAX_POS = 5; # don't consider base-calls close to ends of read 
#MAX_MM = 3;   # maximum number of mismatches allowed for a read to be used, this is computed dynamically using the read length: 36->3, 51->4, 76->5 
#MIN_MAXMQ = 20; # minimum maximum mapping quality of reads required to call SNP
#MIN_READS = 5;  # minimum number of reads covering a position for it to be called as SNP in a sample 
#MIN_PRIOR = 0.001;  # prior of heterozygote, can be changed 
MIN_PRIOR_LOW = 0.2;  # prior of heterozygote once position has been identified as a SNP position 
Qoffset = 33; 

#################################################################################################################################
parser = OptionParser();
parser.add_option("-d","--pileupdir",dest="pileupdir",help="directory with pileup files");
parser.add_option("--rl",dest="readlength",type="int",help="read length, If you have reads of different lengths, e.g. 35 and 40 bp reads, you can specify the longer read length (40) using the --rl option. Alternatively, read length can be specified using the 10th column of the pileup file with the --rl=0 option, default value 36",default=36);
parser.add_option("-c","--maxcycles",type="int",dest="posbins",help="max number of sequencing cycles, this is equal to the readlength for single end reads and 2 x readlength for paired-end reads. For variable length reads, set this to 2 x maximum readlength",default=36);
parser.add_option("--qvoffset",dest="Qoffset",type="int",help="quality value offset, 33 for Sanger, 64 for Illumina base qualities",default=33);
parser.add_option("--mbq",dest="MIN_Q",type="int",help="minimum base quality to consider a base for snp calling, default 10",default=10);
parser.add_option("--mmq",dest="MIN_M",type="int",help="minimum read mapping quality to consider a read for snp calling, default 10",default=10);
parser.add_option("--maxm",dest="MAX_MM",type="int",help="maximum number of mismatches allowed for read to be considered for snp calling. The number of mismatches for each read is stored in column 9 of the pileup file. If this column is missing, set this value to 0. (default value is 3 for 36 bp reads)",default=3);
parser.add_option("--mfile",dest="SNPfile",help="MAQ SNP calls for all samples",default="");
parser.add_option("--dbfile",dest="ncbifile",help="list of dbSNP variants in sequenced region",default="");
parser.add_option("--clipl",dest="MIN_POS",type="int",help="ignore the first x bases of a read for snp calling, default 1",default=1);
parser.add_option("--clipr",dest="MAX_POS",type="int",help="ignore the last x bases of a read for snp calling, default 1",default=1);
parser.add_option("--minmq",dest="MIN_MAXMQ",type="int",help="minimum mapping quality of at least one read covering a variant required to call it as SNP, default 20",default=20);
parser.add_option("--minc",dest="MIN_READS",type="int",help="minimum number of reads required for calling position as SNP, default 5",default=5);
parser.add_option("--hetprior",dest="MIN_PRIOR",type="float",help="heterozygote prior, default 0.001",default=0.001);

(options,args) = parser.parse_args();
pileupfiledir= options.pileupdir; posbins = options.posbins; MIN_Q = options.MIN_Q; MIN_M = options.MIN_M;
MAX_MM = options.MAX_MM; SNPfile = options.SNPfile; ncbifile = options.ncbifile; 
MIN_POS = options.MIN_POS; MAX_POS = options.MAX_POS; MIN_MAXMQ = options.MIN_MAXMQ; MIN_READS = options.MIN_READS; MIN_PRIOR = options.MIN_PRIOR;
readlength = options.readlength; Qoffset = options.Qoffset;

print pileupfiledir,posbins,MIN_Q,MIN_M,MAX_MM,SNPfile,ncbifile,MIN_POS,MAX_POS,MIN_MAXMQ,MIN_READS,MIN_PRIOR,readlength;
print >>sys.stderr,pileupfiledir,posbins,MIN_Q,MIN_M,MAX_MM,SNPfile,ncbifile,MIN_POS,MAX_POS,MIN_MAXMQ,MIN_READS,MIN_PRIOR,readlength;

#sys.exit();


if pileupfiledir == None:
	print '\n-----------------------------------------------------------------------------------------------------';
	print '\n                    SNIP-Seq: program to detect SNPs from Population Sequencing'; 
	print '                      Please enter a path to the directory with pileup files, and the readlength';
	print '\n-----------------------------------------------------------------------------------------------------';
	sys.exit(); 
else:
	filelist = []; filenamelist = []; files =0; 
	if os.path.isdir(pileupfiledir):
		for pileupfile in glob.glob(pileupfiledir + '/*pileup'): 
			filelist.append(open(pileupfile,'r'));  filenamelist.append(pileupfile);  files +=1;
		print >>sys.stderr, 'no of pileup files',files; 
		if files ==0: print >>sys.stderr, 'please provide at least one pileup file'; sys.exit();

	else:	
		for pileupfile in glob.glob(pileupfiledir): filelist.append(open(pileupfile,'r'));  filenamelist.append(pileupfile);  files +=1;
		print >>sys.stderr, 'no of pileup files',files;
		if files ==0: print >>sys.stderr, 'please provide at least one pileup file'; sys.exit();
	print pileupfiledir,SNPfile,ncbifile,files;


#################################################################################################################################

refsequence = ''; positions =0;  postolocus = []; SNPlist_MAQsnps = {};
refsequencearray = array('c');
UniqueTable = []; 

#print >>sys.stderr,'MAQsnpfile',SNPfile,'dbsnpfile',ncbifile;
print >>sys.stderr,'extracting reference sequence from a pileup file';
subprocess.call('awk \'{ print $1,$2,$3 }\' ' + filenamelist[0] + ' > ' + filenamelist[0] + '.3cols',shell=True);

tempfile = open(filenamelist[0]+ '.3cols','r');
for line in tempfile:
	base = line.split();
	refsequencearray.append(base[2]); 
	positions +=1; 	
	if positions%500000==0: print >>sys.stderr, 'done for',positions;
tempfile.close();
subprocess.call('rm -f ' + filenamelist[0] + '.3cols',shell=True);
refsequence += ''.join([refsequencearray[r] for r in xrange(len(refsequencearray))]); 
print >>sys.stderr, 'read in refsequences of total length: ',len(refsequence),'bases';
random.seed();

MAQsnps = snpcalling.read_MAQsnpfile(SNPfile,SNPlist_MAQsnps); 
NCBIlist= {}; snpcalling.read_dbsnpfile(ncbifile,NCBIlist);
if printflag: print >>sys.stderr, 'number of samples',len(filelist),'SNPs',MAQsnps;


##############################################################################################################################

counts = []; indcounts = []; qcounts = []; qsums = []; 	sums = []; mscores = []; 
for r in range(8*posbins+8): counts.append([0,0,0,0]); mscores.append(0);  qcounts.append([0,0,0,0]);
for r in range(2*posbins): qsums.append([0,0,0,0]); # qsums is sum of qcounts for four bases together  
for r in range(2*posbins): sums.append(0.0); # sum of counts 
for j in range(len(filelist)): 
	indcounts.append([]);
	for r in range(8*posbins+8): indcounts[j].append([0,0,0,0]); # individual counts for each bin

readstosample = 255;
genotypematrix = []; # each row is a SNP and each column is an individual 
variantsitelist = []; # first row of this matrix 
noofsnps = 0;

print >>sys.stderr, 'number of samples',len(filelist),'previous SNPs',MAQsnps,'starting to process pileup files for SNPs';
############################################# iterate over all positions in the sequence ##############################################
for i in range(positions):
	if i%100000==0: print >>sys.stderr, 'analyzed pileup file...bases completed',i;
	snpcalling.init_arrays(counts,posbins,indcounts,len(filelist));
	snparray = []; refbase = 0; recalibrated_qarray = []; # recalibrated quality values using population data for this position June 2 09
	coveragelist = [];	potsnplist = []; potsnps = 0;

#	print >>sys.stderr, 'locusname',postolocus[i][0],postolocus[i][1];
	# read in the pileup file for each sample and compute the counts table 
	evaluate =0;
	for j in range(len(filelist)): 
		line = filelist[j].readline(); 
		snp = space.split(line.rstrip('\n')); 
		snparray.append(snp); # list of SNPs lines
		alts = snpcalling.filter_snps(snp);
		if alts >=3: evaluate =1;
	postolocus = [snparray[0][0],snparray[0][1]];
	try:  
		ncbi = NCBIlist[snparray[0][0] + ':' + snparray[0][1]];
		dbsnpinfo = [1,ncbi];
	except KeyError: dbsnpinfo = [0,'-/-'];
	try: MAQsnpinfo = SNPlist_MAQsnps[snparray[0][0] + ':' + snparray[0][1]];
	except KeyError: MAQsnpinfo = 0;

	if evaluate ==0: continue; # not a SNP 

	for j in range(len(filelist)): 
		snp = snparray[j];
		refbase = snpcalling.bti(snp[2]);
		snpinfo = snpcalling.compute_counts(snp,counts,qsums,j,refbase,indcounts,posbins,filenamelist[j],qcounts,0,MIN_Q,MIN_M,MAX_MM,readlength,Qoffset);
		if snpinfo[0] == 1: 
			ref =  refbase; alt = snpinfo[1];
# need to handle SNPs with low mapping quality bases 
			[var,base1,altbase,llf,llr,g11f,g11r,readcounts,snpinfo] = snpcalling.compute_genotypes(snp,filenamelist[j],3,ref,alt,hetratio,3,-1,snp[5],MIN_Q,MIN_M,MIN_POS,MAX_POS,MAX_MM,posbins,readstosample,UniqueTable,i,readlength,Qoffset);

			het0 = 0; het1=0; hetboth = 0; 
			hetprior = math.log(MIN_PRIOR,10); refprior = math.log(0.5-0.5*MIN_PRIOR,10);
			prob00 = refprior; prob01 = llf+llr+hetprior; prob11 = llf+llr+g11f+g11r+refprior;  # three genotype probabilities relative to prob00 
			maxprob = max(prob00,prob01,prob11); 
			prob00 -= maxprob; prob01 -= maxprob; prob11 -= maxprob;
			sumprob = math.log(math.pow(10,prob00) + math.pow(10,prob01) + math.pow(10,prob11),10); 
			prob00 -= sumprob; prob01 -= sumprob; prob11 -= sumprob; 
			sampleprob = random.random(); 
			if sampleprob < math.pow(10,prob00): hetboth = 0;
			elif sampleprob < math.pow(10,prob00)+math.pow(10,prob01): hetboth = 1;
			else: hetboth = 2; 
			potsnplist.append([hetboth,hetboth]); 
			if hetboth > 0: potsnps +=1;
		else: potsnplist.append([0,0]);
		coveragelist.append([0,j]); 
	if potsnps ==0: 
		#if i%100 ==0: print 'no snps',i+1; 
		continue; 
	else:
		print 'initial potential variant site',i+1,potsnps;

#####################################################LOOP OVER TO STabilize###############################################################
	# call function to compute qcounts and qsums using individuals that are not potential variant sites 
	# recompute the new quality scores using qsums and qcounts 
	bt = snpcalling.besttwoQ20bases(counts,posbins); ref = bt[0]; alt = bt[1]; 
	if alt == refbase: tempbase = alt; alt = ref; ref = tempbase; # make refbase to be always the first one
	position_altbase = alt; # global alternate base at this position 
	prev_snps = 0; 
	vars =0; var = 0; hets =0; alts=0; good_snp =0; lastsnps = [0,0]; currsnps = [0,0];
	multiplevars =0; # multiplevars is for two or more individuals with LLR >= 2
	for beta in range(10):
		prev_snps = good_snp; lastsnps = currsnps; 
		for r in range(8*posbins+8): qcounts[r][0] =0;  qcounts[r][1] =0; qcounts[r][2] = 0; qcounts[r][3] = 0;
		for r in range(2*posbins): qsums[r][0] = 0; qsums[r][1] = 0;  qsums[r][2] =0; qsums[r][3] = 0; 
		for r in range(8*posbins):
			for j in range(len(potsnplist)):
				if potsnplist[j][0] ==0 and r%8 < 4: 
					qcounts[r][0] += indcounts[j][r][0]; qcounts[r][1] += indcounts[j][r][1]; qcounts[r][2] += indcounts[j][r][2]; qcounts[r][3] += indcounts[j][r][3];
					qsums[int(r/4)][0] += indcounts[j][r][0]; qsums[int(r/4)][1] += indcounts[j][r][1]; qsums[int(r/4)][2] += indcounts[j][r][2]; qsums[int(r/4)][3] += indcounts[j][r][3];
				if potsnplist[j][0] == 0 and r%8 >= 4: 
					qcounts[r][0] += indcounts[j][r][0]; qcounts[r][1] += indcounts[j][r][1]; qcounts[r][2] += indcounts[j][r][2]; qcounts[r][3] += indcounts[j][r][3];
					qsums[int(r/4)][0] += indcounts[j][r][0]; qsums[int(r/4)][1] += indcounts[j][r][1]; qsums[int(r/4)][2] += indcounts[j][r][2]; qsums[int(r/4)][3] += indcounts[j][r][3];

		recalibrated_qarray = []; good_snp =0; currsnps = [0,0]; snpflag = 0; multiplevars = 0;
		for j in range(len(filelist)): 
			snp = snparray[j]; refbase = snpcalling.bti(snp[2]); 
			MAQSNP = MAQsnpinfo;
			recal_qscores = snpcalling.compute_qvalues(snparray[j],qcounts,counts,qsums,j,posbins,potsnplist,indcounts,refbase,MAQSNP,Qoffset); 
			recalibrated_qarray.append(recal_qscores); pflag = 5; 
			[var1,base11,altbase1,llf1,llr1,g11f1,g11r1,readcounts,snpinfo] = snpcalling.compute_genotypes(snparray[j],filenamelist[j],pflag,ref,alt,hetratio,3,-1,recal_qscores,MIN_Q,MIN_M,MIN_POS,MAX_POS,MAX_MM,posbins,readstosample,UniqueTable,i,readlength,Qoffset);

			het0 = 0; het1=0; hetboth =0;
			hetprior = math.log(MIN_PRIOR,10); 
			refprior = math.log(0.5-0.5*MIN_PRIOR,10);
			prob00 = refprior; prob01 = llf1+llr1+hetprior; prob11 = llf1+llr1+g11f1+g11r1+refprior;  # three genotype probabilities relative to prob00 
			maxprob = max(prob00,prob01,prob11); 
			prob00 -= maxprob; prob01 -= maxprob; prob11 -= maxprob;
			sumprob = math.log(math.pow(10,prob00) + math.pow(10,prob01) + math.pow(10,prob11),10); 
			prob00 -= sumprob; prob01 -= sumprob; prob11 -= sumprob; 
			hetboth =0;
			sampleprob = random.random(); 
			if sampleprob <= math.pow(10,prob00): hetboth = 0;
			elif sampleprob <= math.pow(10,prob00)+math.pow(10,prob01): hetboth = 1;
			else: hetboth = 2; 
			if hetboth > 0: currsnps[0] +=1; currsnps[1] +=1; good_snp +=1; 
			potsnplist[j][0]=hetboth; potsnplist[j][1]=hetboth;
			if hetboth ==1: snpflag = 1; 
			#if hetboth ==1 and prob01 >= 1: multiplevars += 1; 
			if hetboth ==2: snpflag = 1; 

		print 'potential variant site after round',beta,i+1,good_snp;
		if prev_snps == good_snp or (snpflag ==0): break; 

##########################################################################################################################################

	ILL_list = []; MAQgenolist = []; hetlist0 =[]; hetlist1 =[];	hetlist = []; # log-likelihood list and hetlist: if individual is het or alt based on individual reads 
	prev_snps = good_snp;
	vars =0; var = 0; hets =0; alts=0; good_snp =0; altbasegood = '-';
	for j in range(len(filelist)):
		snp = snparray[j]; refbase = snpcalling.bti(snp[2]); 
		[var,base1,altbase,llf,llr,g11f,g11r,readcounts,snpinfo] = snpcalling.compute_genotypes(snp,filenamelist[j],3,ref,alt,hetratio,2,-1,snp[5],MIN_Q,MIN_M,MIN_POS,MAX_POS,MAX_MM,posbins,readstosample,UniqueTable,i,readlength,Qoffset);
		pflag = 1;
		[var1,base11,altbase1,llf1,llr1,g11f1,g11r1,readcounts,snpinfo] = snpcalling.compute_genotypes(snp,filenamelist[j],pflag,ref,alt,hetratio,3,-3,recalibrated_qarray[j],MIN_Q,MIN_M,MIN_POS,MAX_POS,MAX_MM,posbins,readstosample,UniqueTable,i,readlength,Qoffset);
		
		totalreads = 0;
		for g in range(8): totalreads += readcounts[g]; 
		coveragelist[j][0] = totalreads;
		ratiof = float(readcounts[altbase1])/(readcounts[altbase1]+readcounts[refbase]+0.1);
		ratior = float(readcounts[altbase1+4])/(readcounts[altbase1+4]+readcounts[refbase+4]+0.1);
		coveragelist[j].append(llf1);  coveragelist[j].append(llr1); coveragelist[j].append(readcounts); 

		het0 = 0; het1=0; hetboth =0;
		hetprior = math.log(MIN_PRIOR,10); 
		refprior = math.log(0.5-0.5*MIN_PRIOR,10);
		prob00 = refprior; prob01 = llf1+llr1+hetprior; prob11 = llf1+llr1+g11f1+g11r1+refprior;  # three genotype probabilities relative to prob00 
		maxprob = max(prob00,prob01,prob11); 
		if maxprob == prob00: hetboth = 0; 
		elif maxprob == prob01: hetboth = 1; 
		elif maxprob == prob11: hetboth = 2; 

		prob00 -= maxprob; prob01 -= maxprob; prob11 -= maxprob;
		sumprob = math.log(math.pow(10,prob00) + math.pow(10,prob01) + math.pow(10,prob11),10); 
		prob00 -= sumprob; prob01 -= sumprob; prob11 -= sumprob; 
		ILL_list.append([llf1,llr1,llf,llr,base11,altbase1,g11f1,g11r1,ratiof,ratior,snpinfo,totalreads,prob00,prob01,prob11]);
		"""
		if math.pow(10,prob01) >= 0.95*sum: hetboth = 1;
		elif math.pow(10,prob01) + math.prob(10,prob11) >= 0.95*sum: hetboth = 2;
		else: hetboth = 0;
	
		"""
		if hetboth ==1 and (llf1 <= -2 or llr1 <= -2): hetboth = 0; # stranded filter 
		if hetboth > 0:
			if totalreads < MIN_READS: hetboth = 0; # low coverage 
			if (snpinfo[0] ==0 or snpinfo[1] ==0) and (snpinfo[0]+snpinfo[1] <= 4): hetboth = 0; # strand filter 
			if snpinfo[0]+snpinfo[1] < 3: hetboth = 0; # alternate reads filter 
			altreadsmiddle = ILL_list[j][10][6][1] + ILL_list[j][10][7][1]; 
			altreadsends = ILL_list[j][10][6][0] + ILL_list[j][10][7][0] + ILL_list[j][10][6][2] + ILL_list[j][10][7][2];
			if altreadsmiddle < 1 or altreadsmiddle <= 0.1*(altreadsmiddle+altreadsends): hetboth = 0;  # INDEL filter 

		if hetboth ==1: hets +=1; hetlist.append([1,altbase1]); altbasegood = snpcalling.itb(altbase);
		elif hetboth ==2: alts +=1; hetlist.append([2,altbase1]); altbasegood = snpcalling.itb(altbase);
		else: hetlist.append([0,altbase1]);
		if hetboth > 0: vars += 1; good_snp +=1; # good_snp is number of snps...
		MAQgenolist.append([0,0]);
		sample = os.path.basename(filenamelist[j]).split('.')[0] + '.'+ os.path.basename(filenamelist[j]).split('.')[1];
		try: 
			geno = SNPlist_MAQsnps[postolocus[0]+ ':' + postolocus[1] + ':' + sample][1];
			if geno == 'A' or geno == 'C' or geno == 'G' or geno == 'T': MAQgenolist[j][0] = 2;
			else: MAQgenolist[j][0] = 1;
		except KeyError: MAQgenolist[j][0] = 0;

	# genotyping....
	bestPOPLL = [0,0]; bestUR = [0,0]; bestMM = [0,0];
	if good_snp > 0 or multiplevars >= 2:
		print 'SNP is good initially.....',good_snp; 
		good_snp = 0;
		for j in range(len(filelist)):
			het0 = 0; het1=0; hetboth = 0; flag =0;
			if (hetlist[j][1]) != position_altbase and hetlist[j][0] ==0: flag = 1; 
# problematic here for triallelic SNPs 
			hetprior = math.log(MIN_PRIOR_LOW,10); refprior = math.log(0.5-0.5*MIN_PRIOR_LOW,10);
			prob00 = refprior; prob01 = ILL_list[j][0]+ILL_list[j][1]+hetprior; prob11 = ILL_list[j][0]+ILL_list[j][1] + ILL_list[j][6]+ILL_list[j][7]+refprior; # genotype probabilities relative to prob00 
			maxprob = max(prob00,prob01,prob11); 
			#print '%3d %1d %2.1f %2.1f %2.1f %2.1f flag %1d' %(j,hetlist[j][0],maxprob,prob00,prob01,prob11,flag),
			if maxprob == prob00: hetlist[j][0] = 0; 
			elif maxprob == prob01: hetlist[j][0] = 1; 
			elif maxprob == prob11: hetlist[j][0] = 2; 
			if flag ==1: hetlist[j][0] =0;
			if hetlist[j][0] == 1 and prob01 > bestPOPLL[0]+bestPOPLL[1]: bestPOPLL=[ILL_list[j][0],ILL_list[j][1]]; bestUR = [ILL_list[j][10][0],ILL_list[j][10][1]]; bestMM = [ILL_list[j][10][2],ILL_list[j][10][3]];
			if hetlist[j][0] == 2 and prob11 > bestPOPLL[0]+bestPOPLL[1]: bestPOPLL=[ILL_list[j][6],ILL_list[j][7]];
			if hetlist[j][0] > 0: good_snp +=1;
			#print 'FINALCALL',hetlist[j][0];

	if good_snp ==0: print 'SNP did not pass all filters and was failed.....\n'; 

	HWchi_pop = snpcalling.compute_HWE(hetlist);
	HWchi_maq = snpcalling.compute_HWE(MAQgenolist);
	genotypematrix.append([]);
	for j in range(len(filelist)): genotypematrix[noofsnps].append(hetlist[j][0]);
	variantsitelist.append([postolocus[0],postolocus[1],refsequence[i],altbasegood]); 
	noofsnps +=1; 

##################################################################################################################################3
	#if good_snp ==0: if pflag: print >>sys.stderr, 'not detected as SNP by population approach....';
	if (good_snp > 0):
		[maqsnps,bestmaqscore] = snpcalling.bestMAQsnp(postolocus,filenamelist,SNPlist_MAQsnps);
		print 'variant',MAQSNP,
		if good_snp > 0: print '1',
		else: print '0', 
		if dbsnpinfo[0] ==1: print '1',
		else: print '0',
		if printflag: print >>sys.stderr, '%5s %6s %1s %1s popvars %3d maqsnps %3d BEST %3d MCOV %3d bestP %3.1f %3.1f %2d %2d %2d %2d' % (postolocus[0],postolocus[1],snpcalling.itb(refbase),snpcalling.itb(position_altbase),good_snp,maqsnps,bestmaqscore,coveragelist[int(len(coveragelist)/2)][0],bestPOPLL[0],bestPOPLL[1],bestUR[0],bestUR[1],bestMM[0],bestMM[1]),
		if printflag: print >>sys.stderr, 'HWE-pop %3d %3d %3d %2.2f' % (HWchi_pop[0],HWchi_pop[1],HWchi_pop[2],HWchi_pop[3]); #HWchi_maq[0],HWchi_maq[1],HWchi_maq[2],HWchi_maq[3]);
		print '%5s %6s %1s %1s popvars %3d maqsnps %3d BEST %3d MCOV %3d bestP %3.1f %3.1f %2d %2d %2d %2d' % (postolocus[0],postolocus[1],snpcalling.itb(refbase),snpcalling.itb(position_altbase),good_snp,maqsnps,bestmaqscore,coveragelist[int(len(coveragelist)/2)][0],bestPOPLL[0],bestPOPLL[1],bestUR[0],bestUR[1],bestMM[0],bestMM[1]),
		print 'HWE-pop %3d %3d %3d %2.2f' % (HWchi_pop[0],HWchi_pop[1],HWchi_pop[2],HWchi_pop[3]),
		if dbsnpinfo[0] ==1: print 'db 1',dbsnpinfo[1],
		else: print 'db 0 -/-',
		print;

		print 'TABLE',i+1,snpcalling.flankingseq(refsequence,i);	snpcalling.print_bintable_vertical(counts,posbins); print '\n';
		coveragelist.sort();
		for b in range(len(filelist)): 
			j = coveragelist[b][1];
			altreadsends = ILL_list[j][10][6][0] + ILL_list[j][10][7][0] + ILL_list[j][10][6][2] + ILL_list[j][10][7][2];
			altreadsmiddle = ILL_list[j][10][6][1] + ILL_list[j][10][7][1];
			sample = os.path.basename(filenamelist[j]).split('.')[0] + '.'+ os.path.basename(filenamelist[j]).split('.')[1];
			print 'GT %20s %6d %12s' %(postolocus[0],int(postolocus[1]),sample),
			if MAQSNP  > 0: print 'Y',
			else: print 'N',
			if good_snp > 0: print 'Y',
			else: print 'N',
			if hetlist[j][0] == -1: print 'M',MAQgenolist[j][0],'P','-','db',
			else: print 'M',MAQgenolist[j][0],'P',hetlist[j][0],'db',
			print dbsnpinfo[0],dbsnpinfo[1],
			print snpcalling.itb(ILL_list[j][4]),snpcalling.itb(ILL_list[j][5]),'|',
			for r in range(4): print '%2d %2d' % (coveragelist[b][4][r],coveragelist[b][4][r+4]),
			print '%3d | %+3.1f %+3.1f | %+3.1f %+3.1f | %0.2f %0.2f | %0.2f %0.2f %0.2f ' %(coveragelist[b][0],ILL_list[j][0],ILL_list[j][1],ILL_list[j][2],ILL_list[j][3],ILL_list[j][8],ILL_list[j][9],ILL_list[j][12],ILL_list[j][13],ILL_list[j][14]),
			if hetlist[j][0] > 0 or MAQgenolist[j][0] > 0: 
				for k in range(4):		print '%2d ' % (ILL_list[j][10][k]),
				print 'POS',altreadsends,altreadsmiddle,
				if altreadsmiddle < 1 and altreadsends + altreadsmiddle >= 3: print 'INDEL',
				elif altreadsmiddle >= 1 and altreadsends + altreadsmiddle >= 3: print 'SNP',
				else: print 'MISSED',
				#print ILL_list[j][10][6][0] + ILL_list[j][10][7][0] + ILL_list[j][10][6][2] + ILL_list[j][10][7][2],
				#print ILL_list[j][10][6][1] + ILL_list[j][10][7][1],
				#if ILL_list[j][10][4] < MIN_DIST_INDEL and ILL_list[j][10][5] < MIN_DIST_INDEL: print 'INDEL',
				#else: print 'SNP',
			if len(snparray[j]) == 10: print 'idpileup',snparray[j][9], 
			try: 
				msinfo = SNPlist_MAQsnps[postolocus[0]+ ':'+ postolocus[1] + ':' + sample];
				print 'MAQ',msinfo[1],msinfo[2];
			except KeyError: print;

for j in range(len(filelist)): filelist[j].close();

