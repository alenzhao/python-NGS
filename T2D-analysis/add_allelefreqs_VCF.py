#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler

#last modified jan 23 2014

#cat ESP6500SI.allchroms.vcf | awk '{ split($8,F,";"); split(F[5],G,"="); split(G[2],A,","); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\tAF="A[2]";"F[3]; }' > ESP6500SI.allchroms.vcf.AA_AF 
#cat ESP6500SI.allchroms.vcf | awk '{ split($8,F,";"); split(F[5],G,"="); split(G[2],A,","); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\tAF="A[2]/100";"F[3]; }' > ESP6500SI.allchroms.vcf.AA_AF 
#python add_allelefreqs_VCF.py AF ../Variant-call-datasets/2000DanishExomes/2000exomes.AF.vcf phase1-190pools/combined.vcf.annotated

# grep VP= plist.515.ASN.CRISP.012514 | cut -f1-8 | grep -v LowDepth | awk '{ split($8,G,";"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"G[8]";"$8; }' > /projects/stsi/vbansal/T2D-pooledseq-july2012/AlleleFreq-data/plist.515.ASN.CRISP.012514.VCF

## script that takes a file with allele freqs and a VCF file and outputs all variants that overlap the regions in the bed file 

AFtable = {}; vars=0; # list of positions to print, others to omit 

File = open(sys.argv[1]);
for line in File:
        if line[0] == '#': continue;
        var = line.strip().split('\t'); INFO = var[7].split(';'); AF = INFO[0].split('=')[1].split(',')[0];
        AFtable[(var[0],int(var[1]))] = [var[3],var[4],float(AF),1.0-float(AF)]
        vars +=1;
File.close();
print >>sys.stderr, "read",vars,"variants from file",sys.argv[1];


variants=0; overlapping=0;
File = open(sys.argv[2]); # pooled VCF file 
for line in File:
        if line[0] == '#':
                #print line,
                continue;
        var = line.strip().split('\t'); chrom = var[0]; position = int(var[1]); alleles = var[4].split(',');
        if 'exonic' not in var[2] and 'splicing' not in var[2]: continue;
        try:
                varinfo = AFtable[(chrom,position)];
                overlapping +=1;
                if varinfo[0] == var[3] and varinfo[1] == alleles[0]: print varinfo[2],
                elif varinfo[0] == var[3] and len(alleles) > 1 and varinfo[1] == alleles[1]: print varinfo[2],
                elif varinfo[1] == var[3] and varinfo[0] == alleles[0]: print varinfo[3],
                elif varinfo[1] == var[3] and len(alleles) > 1 and varinfo[0] == alleles[1]: print varinfo[3],
                else: print '0.0',

                """
                if varinfo[0] == var[3] and varinfo[1] == alleles[0]: print 'MATCH:'+ varinfo[1]+ ':'+ varinfo[0] + ':' + `varinfo[2]`,
                elif varinfo[0] == var[3] and len(alleles) > 1 and varinfo[1] == alleles[1]: print 'MATCH2:'+ varinfo[1]+ ':'+ varinfo[0] + ':' + `varinfo[2]`,
                elif varinfo[1] == var[3] and varinfo[0] == alleles[0]: print 'MATCH_minor:'+ varinfo[0]+ ':'+ varinfo[1] + ':' + `varinfo[3]`,
                elif varinfo[1] == var[3] and len(alleles) > 1 and varinfo[0] == alleles[1]: print 'MATCH2_minor:'+ varinfo[1]+ ':'+ varinfo[0] + ':' + `varinfo[3]`,
                """
                #else: print 'posmatch:'+varinfo[0] + ':'+ varinfo[1] + ':'+varinfo[2],
        except KeyError: print '0.0',

        print var[0],var[1],var[2],var[3],var[4],var[5],var[6],var[7]
        #print
        variants +=1;
File.close();

print >>sys.stderr, overlapping,'of',variants,"variants overlap with bed file";
~                                                                                                                                                                                                                                                                               
~                                                                                                                                                                                                                                                                               
~                                                                                                                                                                                                                                                                               
~                                                                                                                                                                                                                                                                               
~                                                                                                                                                                                                                                                                               
~                                                                                          
