#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler
import random

#last modified jan 23 2014

## do not handle tri-allelic variants for now... best for such variants would be to split them prior to annotation

## python filter_rarevariants_VCF.py 6pops.allelefreq 190pools.ancestry.5pops combined.vcf.annotated > b

if len(sys.argv) < 4: print >>sys.stderr, "python filter_rarevariants_VCF.py 6pops.allelefreq 190pools.ancestry.5pops combined.vcf.annotated 0/1(printVCF)"; sys.exit();

random.seed();

AFtable = {}; vars=0; # list of positions to print, others to omit 
printVCF = 0;
if len(sys.argv) > 4: printVCF = int(sys.argv[4])

# EUR.1KG.505 ASN.1KG.515 AFR.1KG.507 EA.4400 AA.2250 Danish.2000
File = open(sys.argv[1]); ## allele frequency file
for line in File:
        var = line.strip().split(); firstcol = 6;
        if line[0] == '#':
                #print var; 
                for i in xrange(len(var)):
                        if var[i] == 'chr' or var[i] == 'chrom': firstcol = i;
                header = var;
                continue;
        else:
                AFtable[(var[firstcol],var[firstcol+1],var[firstcol+3],var[firstcol+4])] = [float(x) for x in var[0:firstcol]]
                vars +=1;
File.close();
print >>sys.stderr, "read",vars,"variants from file",sys.argv[1];

## read ancestry file for the pools 

#pool AFR   ASN   EUR   MEX Gujarati
ancestryINFO = {};
File = open(sys.argv[2]); ## 
for line in File:
        if line[0] == '#': continue;
        pool = line.strip().split();
        ancestryINFO[pool[0]] = [float(x) for x in pool[1:]]
File.close();

## if it is difficult to decide, then set genotype to missing, otherwise set to 0, and reduce effective pool size for that variant only...

variants=0; overlapping=0; poollist = [];
File = open(sys.argv[3]); # pooled VCF file 
for line in File:
        if line[0] == '#' and line[1] == 'C' and line[2] == 'H': var = line.strip().split('\t'); poollist = var;

        if line[0] == '#':
                if printVCF ==1: print line,
                continue;

        var = line.strip().split('\t'); chrom = var[0]; position = int(var[1]); alleles = var[4].split(',');
        if 'exonic' not in var[2] and 'splicing' not in var[2]: continue;
        if 'nonsyn' not in var[2] and 'frameshift' not in var[2] and 'splicing' not in var[2] and 'stop' not in var[2]: continue;
        try:
                varinfo = AFtable[(var[0],var[1],var[3],var[4])]; #print varinfo;
                potential = 0;
                if varinfo[1] > varinfo[0] or varinfo[2] > varinfo[0] or varinfo[4] > varinfo[3]: potential = 1;
                if varinfo[0] > 0.01: potential = 0;

                if varinfo[3] < 0.000001: varinfo[3] = 0.00005;
                if varinfo[5] < 0.000001: varinfo[5] = 0.000125;

                if varinfo[0] < 0.000001: varinfo[0] = varinfo[5];
                if varinfo[0] < 0.000001: varinfo[0] = 0.0005;


                overlapping +=1;
                if len(alleles) > 1 or potential ==0:
                        if printVCF==1: print line,
                else:
                        if printVCF >= 1:
                                for i in xrange(8): sys.stdout.write(var[i] + '\t')
                                sys.stdout.write(var[8]);

                        #if position== 17408485: print 'potential',potential,position;
                        for i in xrange(9,len(var)):
                                origin = 'EUR';
                                PS = [varinfo[0]*ancestryINFO[poollist[i]][2],varinfo[1]*ancestryINFO[poollist[i]][1],varinfo[2]*ancestryINFO[poollist[i]][0]];
                                sum = PS[0] + PS[1] + PS[2];
                                if sum > 0: PS[0] /= sum; PS[1] /= sum; PS[2] /= sum;
                                PS1 = [varinfo[3]*ancestryINFO[poollist[i]][2],varinfo[5]*ancestryINFO[poollist[i]][2],varinfo[4]*ancestryINFO[poollist[i]][0]];
                                sum = PS1[0] + PS1[1] + PS1[2];
                                if sum > 0: PS1[0] /= sum; PS1[1] /= sum; PS1[2] /= sum;

                                if varinfo[1] > varinfo[0] and ancestryINFO[poollist[i]][1] > 0: origin = 'ASN';
                                elif varinfo[2] > varinfo[0] and ancestryINFO[poollist[i]][0] > 0: origin = 'AFR1';
                                elif varinfo[4] > varinfo[3]  and varinfo[4] > varinfo[5] and ancestryINFO[poollist[i]][0] > 0: origin = 'AFR2';
                                genotype = var[i].split(':');
                                filter_genotype = 0; gt = int(genotype[0]);
                                threshold = 0.5;

                                if origin != 'EUR' and gt > 0 and gt <=2 and PS[1] >= threshold: filter_genotype = 1;
                                elif origin != 'EUR' and gt > 0 and gt <=2 and PS[2] >= threshold and PS1[2] >=threshold-0.1: filter_genotype = 3;
                                elif origin != 'EUR' and gt >0 and gt <=2 and PS[2] >= threshold-0.1 and PS1[2] >= threshold: filter_genotype = 4;
                                elif origin != 'EUR' and gt > 0 and gt <=2 and PS[1]+PS[2] >= 0.9: filter_genotype = 2;

                                """     
                                if origin != 'EUR' and int(genotype[0]) > 0 and random.random() < PS[1]: filter_genotype = 1; ## Asian allele 
                                elif origin != 'EUR' and int(genotype[0]) > 0 and 2*random.random() < PS[2]+PS1[2]: filter_genotype = 1;  ## african alleles
                                """

                                if printVCF >= 1 and filter_genotype ==0: print "\t%s" %(var[i]),
                                elif printVCF >= 1 and filter_genotype >=1:
                                        if genotype[0] == '1': newvar = ':'.join(['0'] + genotype[1:] + ['FILTER']);
                                        else: genotype[0] = '.'; newvar = ':'.join(['.'] + genotype[1:] + ['FILTER']);
                                        #genotype[0] = 'FILTER'; newvar = ':'.join(genotype);
                                        sys.stdout.write("\t" + newvar);

                                if (printVCF ==0 and origin != 'EUR' and int(genotype[0]) > 0 and filter_genotype >=1):# or (position == 17408485 and int(genotype[0]) > 0): 
                                        print var[0],var[1],var[3],var[4],var[i],varinfo,poollist[i],ancestryINFO[poollist[i]],origin,var[2];
                                        print 'EUR:',PS[0],'ASN:',PS[1],'AFR:',PS[2],'|','EUR1:',PS1[0],'EUR2:',PS1[1],'AA:',PS1[2];
                                        print 'fil',filter_genotype;

                        if printVCF >= 1: print '\n',
        except KeyError:
                if printVCF ==1: print line,
                print >>sys.stderr,  var[0],var[1],var[2],var[3],var[4]
        variants +=1;

File.close();

print >>sys.stderr, overlapping,'of',variants,"variants overlap with bed file";
