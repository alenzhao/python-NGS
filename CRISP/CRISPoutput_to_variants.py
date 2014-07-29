#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math


#### author: VIKAS BANSAL Sept 2 2010 #####

if len(sys.argv) < 4: print 'input arguments: CRISP.out number.samples haplotypes.per.pool'; sys.exit();

File = open(sys.argv[1]); samples = int(sys.argv[2]); ploidy = int(sys.argv[3]);

#############################################################################################################

print '##format=VCFlike';
print '##source=CRISP.variantcaller';
print '##QUAL=CTpvalue;QVpvalue_fwd;QVpvalue_rev';
print '##NS=number of samples';
print '##VS=samples with alternate allele';
print '##GT=Genotype';
print '##GQ=Genotype Quality, Phred scaled consensus score';
print '##DP=Depth of Coverage';
print '##A1=allele1';
print '##A2=allele2';
print '##A1C=reads.with.allele1';
print '##A2C=reads.with.allele2';
print "%s\t%9s\t%6s\t%1s\t%1s\t%20s\t%5s\t%20s\t%10s\t" %('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'),
#############################################################################################################


### handling multi-allelic variants, genotypes reported as 1/3 where 1 and 3 refer to the two alleles ######

#### for non-diploid samples, we want to report the allele counts based on pool size ##### 

#### alternate format, report information for each genotype together #####

###### TO FIX: false indels (third allele), homopolymer indels using MCMC algorithm, reference allele bias in favor of indels....


snps =0; indels =0; variants =0;

while 1:
	line = File.readline();
	if not line: break;
	var = line.strip().split(); l = len(var); vars =0;
	for i in xrange(l):
		if var[i] == 'variants': vars = int(var[i+1]); break; 

	if vars > 0:
		type = var[i+2]; first = i+5; variants +=1;
		CTpval = var[15]; QVpvalf = var[19]; QVpvalr = var[20]; poolf = var[21]; poolr = var[22];
		coverage = int(var[8]) + int(var[10]); 
		allele0 = var[4]; allele1 = ''; allele2 = ''; type = 'snp';
		if var[5] == var[4]: 
			allele1 = var[6]; altallele = var[6];
			if '+' in var[6] or '-' in var[6]: type = 'indel';
		elif var[6] == var[4]: 
			allele1 = var[5]; altallele = var[5]; 
			if '+' in var[5] or '-' in var[5]: type = 'indel';
		else: 
			allele1 = var[5]; allele2 = var[6];  altallele = var[5] + ',' + var[6];
			if '+' in var[5] or '-' in var[5]: type = 'indel';

		if type == 'snp': snps +=1; number = snps;
		if type == 'indel': indels +=1; number = indels;
		varinfo = type + ':' + `number`;

		genotypes = []; phredscores = []; allelecounts = []; samplelist = [];
		for s in xrange(samples):
			gen = File.readline().split(); allelecounts.append([gen[6],gen[7]]); 
			genotypes.append(gen[8]); phredscores.append(gen[9]);
			if variants ==1: samplelist.append(gen[11].split('/')[-1].strip('.pileup'));

		if variants ==1:
			for s in xrange(samples): print samplelist[s],
			print;

		#if var[0] != 'non': varinfo += ';' + var[0] + ':' + var[1];
		#else: varinfo += ';novel';
		if type == 'indel' and float(CTpval) >= -2: # and float(QVpvalf)  >= -5 and float(QVpvalr) >= -5 : 
			print '%s\t%9s\t%6s\t%1s\t%1s\t%30s\t%5s\t%20s\t' %(var[2],var[3],varinfo,var[4],altallele,CTpval+';'+QVpvalf+';'+QVpvalr+ ';' + poolf + ';' + poolr,'HP_indel','NS='+`samples`+';DP='+`coverage`+ ';VS=' + `vars`),
		else:
			print '%s\t%9s\t%6s\t%1s\t%1s\t%30s\t%5s\t%20s\t' %(var[2],var[3],varinfo,var[4],altallele,CTpval+';'+QVpvalf+';'+QVpvalr + ';' + poolf + ';' + poolr,'PASS','NS='+`samples`+';DP='+`coverage`+ ';VS=' + `vars`),
			

		if ploidy <= 2:
			print 'GT:GQ:DP\t',
			for s in xrange(samples):
				if var[5] == var[4]: 
					if genotypes[s] == '0': gt = '0/0';
					if genotypes[s] == '2': gt = '1/1';
					if genotypes[s] == '1': gt = '0/1';
				elif var[6] == var[4]: 
					if genotypes[s] == '0': gt = '1/1';
					if genotypes[s] == '2': gt = '0/0';
					if genotypes[s] == '1': gt = '0/1';
				elif var[6] != var[4]: 
					if genotypes[s] == '0': gt = '2/2';
					if genotypes[s] == '2': gt = '1/1';
					if genotypes[s] == '1': gt = '1/2';
					
				print gt + ':' + phredscores[s] + ':' +  `int(allelecounts[s][0]) + int(allelecounts[s][1])`,
				#print genotypes[s] + ':' + phredscores[s] + ':' +  allelecounts[s][0] + ':' + allelecounts[s][1],
			print;
		else:
			print 'A1C:A2C\t',
			print ' '.join([allelecounts[s][0] + ':' + allelecounts[s][1] for s in xrange(samples)]);
			
			"""
			print ''.join([genotypes[s] for s in xrange(samples)]),
			print ','.join([phredscores[s] for s in xrange(samples)]),
			print ''.join([genotypes[s] for s in xrange(samples)]),
			print ','.join([allelecounts[s][0] + ':' + allelecounts[s][1] for s in xrange(samples)]);
			"""
			
File.close();	


  

