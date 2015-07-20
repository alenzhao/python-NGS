#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math


snps =0; indels =0; variants =0;

print '##format=VCFlike';
print '##source=CRISP.pooled.variantcaller';
print '##QUAL=CTpvalue;QVpvalue_fwd;QVpvalue_rev';
print '##NP=number.of.pools';
print '##DP=total.depth.across.all.pools';
print '##VP=pools.with.alternate.allele';
print '##A1=allele1';
print '##A2=allele2';
print '##A1C=reads.with.allele1';
print '##A2C=reads.with.allele2';
print "%20s\t%9s\t%6s\t%1s\t%1s\t%20s\t%5s\t%20s\t%10s\t" %('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'),

File = open(sys.argv[1]);
for line in File:
	var = line.strip().split(); l = len(var); vars =0;
	for i in xrange(l):
		if var[i] == 'variants': vars = int(var[i+1]); break; 

	if vars > 0:
		type = var[i+2]; first = i+5; samples  = l-first;
		if variants ==0: # print first line 
			for i in xrange(first,l): s = var[i].split(':'); print s[0]+'\t',
			print;
		
		variants +=1;
		if type == 'snp': snps +=1; number = snps;
		if type == 'indel': indels +=1; number = indels;
		CTpval = var[15]; QVpvalf = var[19]; QVpvalr = var[20];
		coverage = int(var[8]) + int(var[10]); 
		if var[5] == var[4]: altallele = var[6];
		elif var[6] == var[4]: altallele = var[5];
		else: altallele = var[5] + ',' + var[6];
		varinfo = type + ':' + `number`;
		if var[0] != 'non': varinfo += ';' + var[0] + ':' + var[1];
		else: varinfo += ';novel';

		print '%20s\t%9s\t%6s\t%1s\t%1s\t%20s\t%5s\t%20s\t%10s\t' %(var[2],var[3],varinfo,var[4],altallele,CTpval+';'+QVpvalf+';'+QVpvalr,'0','NP='+`samples`+';DP='+`coverage`+ ';VP=' + `vars` + ';A1='+var[5] + ';A2='+var[6],'A1C:A2C'),
		for i in xrange(first,l): s = var[i].split(':'); print s[1]+':'+s[2],
		print;

File.close();	


  

