#!/usr/bin/python
import os, glob,sys;
from readfasta import read_fasta_new

### VCF format for indel file, FORMAT = 'DELINS' or FORMAT = 'VCF';
### python leftjustify-indels.py indel-file reference.fasta FORMAT 

def read_indels(indelfile,seqfile,FORMAT):
        fastaindex = {}; sequences = read_fasta_new(seqfile,fastaindex);

        File = open(indelfile,'r');
        for line in File:
                if line[0] == '#': print line.strip(); continue;

                indel = line.split(); chrom = indel[0]; position = int(indel[1]);
                ref = indel[2]; alt = indel[3];
                if length(ref) == length(alt): print line.strip(); continue;
                elif length(ref) < length(alt): type = 'I'; position +=1; ibases = alt[1:];
                elif length(ref) > length(alt): type = 'D'; position +=1; ibases = ref[1:];

                if variant[0] == '-': type = 'I'; ibases = variant[1].upper();
                else: type = 'D'; ibases = variant[0].upper();
                lengthvar = len(ibases);

                index = fastaindex[chrom];
                if type == 'D':
                        leftpos = position-2; rightpos = position+lengthvar-2; leftshift = 0;
                        while sequences[index][1][leftpos] == sequences[index][1][rightpos]:
                                leftpos -=1; rightpos -=1; leftshift +=1;
#                       print chrom,position,type,bases, sequences[index][1][position-1-leftshift:position-1],sequences[index][1][position-1:position-1+lengthvar],leftshift;
                #       print line,
                        if FORMAT == 'DELINS': print 'DEL',chrom,position-leftshift,sequences[index][1][position-1-leftshift:position-1-leftshift+lengthvar];
                        if FORMAT == 'VCF': print chrom,position-leftshift-1,sequences[index][1][position-1-leftshift-1:position-1-leftshift+lengthvar],sequences[index][1][position-1-leftshift-1];
                elif type == 'I':
                        offset = 1000;
                        if position-1-offset < 0: offset = position-1;

                        newseq = sequences[index][1][position-1-offset:position-1] + bases; l = len(newseq);
                        rightpos = l-1; leftpos = l-1-lengthvar; leftshift = 0;

                        while newseq[leftpos] == newseq[rightpos] and leftpos >0: leftpos -=1; rightpos -=1; leftshift +=1;
                        #print chrom,position,type,bases, sequences[index][1][position-1-leftshift:position-1],bases,leftshift;
                #       print line,
                        if FORMAT == 'DELINS': print 'INS',chrom,position-leftshift,newseq[rightpos-lengthvar+1:rightpos+1];
                        if FORMAT == 'VCF': print chrom,position-leftshift-1,newseq[rightpos-lengthvar],newseq[rightpos-lengthvar:rightpos+1];

                continue;


