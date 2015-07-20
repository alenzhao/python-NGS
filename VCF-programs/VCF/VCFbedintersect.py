#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler

#last modified feb 27 2012
## script that takes a BEDfile and a VCF file and outputs all variants that overlap the regions in the bed file 
# python /home/vbansal-scripps/PROGRAMS/INDELS-Exome-scripts/extract_vcf_bedfile.py /home/vbansal-scripps/GENOMES/Nimblegen.exome/Nimblegen.exome.chr20.50flank.bed BGIdata/LuCAMP_200exomeFinal.maf.chr20  > BGIdata/LuCAMP_200exomeFinal.maf.chr20.50bpoverlap

def VCFintersect(vcffile,bedfile):
        pileuptable = {}; # list of positions to print, others to omit 
        File = open(bedfile,'r');
        for line in File:
                var = line.split();
                for pos in xrange(int(var[1]),int(var[2])+1): pileuptable[(var[0],pos)] = 1;
        File.close();
        #print >>sys.stderr, "read",indels,"indels from file",sys.argv[1];

        variants=0; overlapping=0;
        File = open(vcffile,'r');
        for line in File:
                if line[0] == '#':
                        print line,
                        continue;
                var = line.split();chrom = var[0]; position = int(var[1]);
                try:
                        bedpresent = pileuptable[(chrom,position)];
                        variants +=1; overlapping +=1;
                        print line,
                except KeyError: variants +=1;
        File.close();
        print >>sys.stderr, overlapping,'of',variants,"variants overlap with bed file";
