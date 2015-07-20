#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler
import random

#last modified aug 7 2014

## do not handle tri-allelic variants for now... best for such variants would be to split them prior to annotation

## python filter_rarevariants_VCF.py 6pops.allelefreq 190pools.ancestry.5pops combined.vcf.annotated > b

#if len(sys.argv) < 4: print >>sys.stderr, "python filter_rarevariants_VCF.py 6pops.allelefreq 190pools.ancestry.5pops combined.vcf.annotated 0/1(printVCF)"; sys.exit();

## arg1 = phase1_pools, arg2 = phase2_pools, arg3 = phase3_pools 
## arg4 = VCF_phase3_44, arg

# python find_carriers.py PoolInfo/phase1_cases PoolInfo/phase2_cases PoolInfo/phase3_pools.77 phase1.CRISP.vcf.annotated  phase2.CRISP.vcf.annotated plist.78pools.VCF.080814 > b

# python find_carriers.py PoolInfo/phase1_cases PoolInfo/phase2_cases PoolInfo/phase3_pools phase1.CRISP.vcf.annotated  phase2.CRISP.vcf.annotated plist.44pools.VCF.080614 > b

##############################################################################################################################

def make_pool_table(samplefile):
	
	PHASE1 = {};
	File = open(samplefile); # phase1_pools 
	for line in File: 
		if line[0] == '#': continue
		pool = line.strip().split();
		for i in xrange(1,len(pool)): PHASE1[pool[i]] = pool[0];
		PHASE1[pool[0]] = [0,pool[1:],[]];
		#print pool[0],pool[1:]
	File.close(); 
	return PHASE1


## read VCF files for phase1 and phase2, all pools have same size...
## need to handle tri-allelic variants 
def readVCF(VCFfile):

	samples = {}; 	varlist = {};  variants =0; samplelist = [];
	File = open(VCFfile,'r');
	for line in File: 
		if line[0] == '#' and line[1] == '#': continue
		elif line[0] == '#' and line[1] == 'C': 
			s = line.strip().split(); 
			for i in xrange(9,len(s)): samples[s[i]] = i; 
			samplelist = s;
		else: 
			var = line.strip().split();
			chrom = var[0]; position = int(var[1]); ref = var[3]; alt = var[4].split(','); 
			for allele in alt: varlist[(chrom,position,allele)]  = var; 
			variants +=1; 
	File.close();
	print >>sys.stderr, 'variants',variants,'in VCF file',VCFfile;
	return [samples,varlist,samplelist] 
			
##########################################################################################


PHASE1 = make_pool_table(sys.argv[1]); ## phase 1 case pools
PHASE2 = make_pool_table(sys.argv[2]); ## phase 1 case pools
#for p in POOL_samples.iterkeys(): print p, POOL_samples[p]; sys.exit();

#print PHASE1['G1P2']

missing = 0;
PHASE3 = {}; PHASE3_1 = {};
PAIRWISE_TABLE = {}; 
File = open(sys.argv[3]); # phase1_pools 
for line in File: 
	if line[0] == '#': continue
	pool = line.strip().split();
	#PHASE3[pool[0]] = pool[1:]; 
	#print pool;
	PHASE3[pool[0]] = [];
	aoo = 0; noo = 0; obese = 0; slist = []; minA= 100; maxA = 0; phase1 = 0;
	#print pool[0],
	for i in xrange(1,len(pool)): 
		if pool[i] in PHASE1: pold = PHASE1[pool[i]]; PHASE3[pool[0]].append([pold,pool[i]]); PHASE1[pold][0] +=1; PHASE1[pold][2].append(pool[i]);
		elif pool[i] in PHASE2: pold = PHASE2[pool[i]]; PHASE3[pool[0]].append([pold,pool[i]]); PHASE2[pold][0] +=1; PHASE2[pold][2].append(pool[i]);
		else: pold = 'missing'; missing +=1;
		if 'obes' in pold: obese +=1; 
		sid = pold.split('(')[0].split(':');
		if len (sid[0].split('_')) == 2 and int(sid[0].split('_')[1]) < 4200: phase1 +=1;
		slist.append(pold); 
		if len(sid) ==2:
			AGE = int(sid[1].split('=')[1]);
			aoo += AGE; noo += 1;
			if AGE < minA: minA = AGE; 
			if AGE > maxA: maxA = AGE; 
		#print sid,
	#print pool[0],len(pool)-1,`minA` + '-' + `maxA`,aoo/noo,phase1,'\t'.join([slist[i] for i in xrange(len(slist))]);
	#print len(pool)-1-phase1
	#print '\n',
File.close();	
#print >>sys.stderr, 'missing',missing;
#sys.exit()

#for pools in PHASE3: print pools,PHASE3[pools];


########################################################################################

ULM_ids = {}; ## 54_H_1 -> Scripps_Ulm_0544
File = open('PoolInfo/allsamples_UlmIDs_withAO'); 
for line in File: s = line.strip().split(); ULM_ids[s[1]] = s[0]; 
File.close(); 

extra_genes = {}; ## list of genes sequenced in phase2/3 but not in phase1, any variant in such gene cannot always be carrier deconvoluted
File = open('PoolInfo/phase2_extragenes'); 
for line in File: g = line.strip().split()[0]; extra_genes[g] = 1; 
File.close();
#for g in extra_genes.iterkeys(): print g;
##############################################################################################################################

#class VCF: def __init__(self):
	

phase1_VCF = readVCF(sys.argv[4]); ## VCF from 190 pools, phase1
phase2_VCF = readVCF(sys.argv[5]); ## VCF from 132 pools, phase2

############### read VCF files for phase1 and phase2, all pools have same size...###############




def print_var_pools(var,samplelist,PHASE): ## VCF file
	for i in xrange(9,len(var)): 
		G = var[i].split(':'); v=0;
		if G[0] != '.' and ',' not in G[0]: v = int(G[0]); #altcount += v; 
		if v > 0 and samplelist[i] in PHASE: print var[i],samplelist[i],PHASE[samplelist[i]];
	print;

def get_var_pools(var,samplelist,PHASE,allele):
	l = []; lf = []; missing = 0;

	alleles = var[4].split(',');
	if alleles[0] == allele: a = 0; 
	elif len(alleles) > 1 and alleles[1] == allele: a =1; 
	elif len(alleles) > 2 and alleles[2] == allele: a =2; 
	else: a = -1; 

	for i in xrange(9,len(var)): 
		G = var[i].split(':'); v=0;
		if G[0] != '.' and ',' not in G[0] and a >=0: v = int(G[0]); #altcount += v; 
		elif G[0] != '.' and ',' in G[0] and a >=0: v = int(G[0].split(',')[a]); #altcount += v; 

		if v > 0 and samplelist[i] in PHASE: l.append(samplelist[i]); 
		elif v >0: missing +=1;
		if v > 0 and samplelist[i] in PHASE: lf.append([samplelist[i],var[i]]); 
	return [l,lf,missing];



## open VCF file for phase3 pools ###################

File = open(sys.argv[6]); # phase3 VCF file 
samplelist = []; 	varlist = {};  variants =0;
for line in File: 
	if line[0] == '#' and line[1] == '#': continue
	elif line[0] == '#' and line[1] == 'C': 
		samplelist = line.strip().split();
	else: 
		var = line.strip().split();
		chrom = var[0]; position = int(var[1]); ref = var[3]; alt = var[4].split(','); 
		if len(alt) > 2: continue; ## more than tri-allelic

		for a in xrange(len(alt)): ## process each allele separately to identify carriers
			allele = alt[a]; 
			altcount =0; refcount =0;
			for i in xrange(9,len(var)): 
				G = var[i].split(':'); v=0;
				if G[0] != '.' and ',' not in G[0]: ## special case where all pools have same size
					v = int(G[0]);
					if v > 0: altcount +=1;
				elif G[0] != '.' and ',' in G[0] and len(alt) ==1: 
					v = G[0].split(','); 
					if int(v[1]) > 0: altcount +=1; #print G[0],v;
				elif G[0] != '.' and ',' in G[0] and len(alt) ==2: 
					v = G[0].split(','); 
					if int(v[a+1]) > 0: altcount +=1; 

			if altcount <=12 and 'PASS' in var[6]: 

				flag = 0;
				if (chrom,position,allele) in phase1_VCF[1] and 'exonic' in phase1_VCF[1][(chrom,position,allele)][2]: flag += 1; var1 = phase1_VCF[1][(chrom,position,allele)]; samplelist1 = phase1_VCF[2]; 
				if (chrom,position,allele) in phase2_VCF[1] and 'exonic' in phase2_VCF[1][(chrom,position,allele)][2]: flag += 2; var2 = phase2_VCF[1][(chrom,position,allele)]; samplelist2 = phase2_VCF[2];

				if (chrom,position,allele) not in phase1_VCF[1] and (chrom,position,allele) not in phase2_VCF[1]:
					if 'splicing'  in var[2] or 'frame' in var[2] or 'stop' in var[2] or 'nonsyn' in var[2]: 	
						gene = var[2].split(';')[1].split(',')[0]; 
						if gene not in extra_genes:
							print '##variant not found in phase1/2', 
							for i in xrange(8): print var[i],
							print '\n';		

				if flag ==0:	continue; ## variant not observed in previous phase1 or phase2 
			
				printflag =0;	
				l1 = []; l2 = [];
				
				if flag ==1 or flag ==3:
					if 'ncRNA_exonic' in phase1_VCF[1][(chrom,position,allele)][2]: continue; 
					if 'splicing' not in phase1_VCF[1][(chrom,position,allele)][2] and 'frame' not in phase1_VCF[1][(chrom,position,allele)][2] and 'stop' not in phase1_VCF[1][(chrom,position,allele)][2] and 'nonsyn' not in phase1_VCF[1][(chrom,position,allele)][2]: continue; 
					print chrom, position, ref,alt,allele,var[6],var[7],altcount,
					if len(alt) > 1: print 'triallelic';
					else: print ''; 
					printflag = 1;
					[l1,l1f,m1] = get_var_pools(var1,samplelist1,PHASE1,allele);
					fields = var1[7].split(';'); 
					genename = var1[2].split(';')[1].split(',')[0]; annotation = var1[2]
					for f in fields: 
						if f[0] == 'A' and f[1] == 'C': AC1 = f
					print 'var_in_phase1:',var1[4],var1[7].split(';')[1],var1[2],AC1,m1,#var1,samplelist1;
					for pool in l1f: 
						print pool[0]+'|'+pool[1]+'|'+`PHASE1[pool[0]][0]`,
						#if PHASE1[pool][0] < 4: print PHASE1[pool][2],
					print;
					#print_var_pools(var1,samplelist1,PHASE1); 
				
				if flag >=2:
					if 'ncRNA_exonic' in phase2_VCF[1][(chrom,position,allele)][2]: continue; 
					if 'splicing' not in phase2_VCF[1][(chrom,position,allele)][2] and 'frame' not in phase2_VCF[1][(chrom,position,allele)][2] and 'stop' not in phase2_VCF[1][(chrom,position,allele)][2] and 'nonsyn' not in phase2_VCF[1][(chrom,position,allele)][2]: continue; 
					if printflag ==0:
						print chrom, position, ref,alt,allele,var[6],var[7],altcount,
						if len(alt) > 1: print 'triallelic';
						else: print ''; 

					[l2,l2f,m2] = get_var_pools(var2,samplelist2,PHASE2,allele)
					fields = var2[7].split(';'); 
					genename = var2[2].split(';')[1].split(',')[0]; annotation = var2[2]
					for f in fields: 
						if f[0] == 'A' and f[1] == 'C': AC2 = f
					print 'var_in_phase2:',var2[4],var2[7].split(';')[1],var2[2],AC2,m2,#var2,samplelist2;
					for pool in l2f: 
						print pool[0]+'|'+pool[1]+'|'+`PHASE2[pool[0]][0]`,
						#if PHASE2[pool][0] < 4: print PHASE2[pool][2],
					print;
					#print_var_pools(var2,samplelist2,PHASE2); 
				
				for i in xrange(9,len(var)): 
					G = var[i].split(':'); v=0;
					if G[0] != '.' and ',' not in G[0]: v = int(G[0]); #altcount += v; 
					elif G[0] != '.' and ',' in G[0] and len(alt) ==1: v = int(G[0].split(',')[1]);  ## VCF with variable pool sizes
					elif G[0] != '.' and ',' in G[0] and len(alt) ==2: v = int(G[0].split(',')[a+1]);  ## VCF with variable pool sizes
					if v > 0 and samplelist[i] != 'G6': 
						pool = samplelist[i]; matches=0;
						print var[0]+':'+var[1],'POOL',pool,var[i],v,
						for pair in PHASE3[pool]: 
							if pair[0] in l1: matches +=1; print pair[0]+':'+pair[1]+':'+`PHASE1[pair[0]][0]`,		
							elif pair[0] in l2: matches +=1; print pair[0]+':'+pair[1]+':'+`PHASE2[pair[0]][0]`,
						print 'matches:',matches,
						if genename in extra_genes: print 'PHASE2_ONLY_'+':'+genename;
						elif matches ==0: print 'MISSED';
						else: print '\n',

				## print information for output table ###

				if len(alt) >=2: print 'CARRIER_TRI',var[0],var[1],var[3],allele,annotation,
				else: print 'CARRIER',var[0],var[1],var[3],allele,annotation,
				if flag ==1 or flag==3: print AC1,
				else: print '-',
				if flag >=2: print AC2,
				else: print '-',
				fields = var[7].split(';');  
				for f in fields: 
					if f[0] == 'A' and f[1] == 'C': AC3 = f
				print AC3,
				
				for i in xrange(9,len(var)): 
					G = var[i].split(':'); v=0;
					if G[0] != '.' and ',' not in G[0]: v = int(G[0]); #altcount += v; 
					elif G[0] != '.' and ',' in G[0] and len(alt) ==1: v = int(G[0].split(',')[1]);  ## VCF with variable pool sizes
					elif G[0] != '.' and ',' in G[0] and len(alt) ==2: v = int(G[0].split(',')[a+1]);  ## VCF with variable pool sizes
					if v > 0 and samplelist[i] != 'G6': 
						pool = samplelist[i]; mlist = [];
						for pair in PHASE3[pool]: 
							if pair[0] in l1: mlist.append(ULM_ids[pair[1]]);
							elif pair[0] in l2: mlist.append(ULM_ids[pair[1]]); 
						if len(mlist) ==0 and genename in extra_genes: print 'no_data',
						elif len(mlist) ==0 : print 'missed',
						elif len(mlist) ==1: print mlist[0],
						elif len(mlist) ==2 and v ==2: print mlist[0],mlist[1],
						elif len(mlist) ==2 and v ==1 and mlist[0] == mlist[1]: print mlist[0],
						elif len(mlist) ==3 and v ==3: print mlist[0],mlist[1],mlist[2],
						elif len(mlist) ==4 and v ==4: print mlist[0],mlist[1],mlist[2],mlist[3],
						elif len(mlist)  > v: print '/'.join([mlist[i] for i in xrange(len(mlist))]), 
			
				print;	
				

				#print_var_pools(var,samplelist,PHASE3); 
				print '---------------------------------------------------\n';	



