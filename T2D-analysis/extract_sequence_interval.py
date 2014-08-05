#!/usr/bin/python2.4
import sys, os, glob, subprocess, time, math, random
from readfasta import read_fasta
###################################################################################################################################

def extract_sequences(fastafile,targetinterval):
        [sequences,reflist] = read_fasta(fastafile);

        revcomp = {'A':'T', 'T':'A','C':'G','G':'C','N':'N','n':'N'};

        chrom = targetinterval.split(':')[0]; start = targetinterval.split(':')[1].split('-')[0];
        end = targetinterval.split(':')[1].split('-')[1];
        targetlist = [[chrom,int(start.replace(',','')),int(end.replace(',',''))]];
        targets =len(targetlist);

        # need to take care of 0/1 offset, 1 corresponds to 0 in sequence....
        for i in xrange(targets):
                chrom = targetlist[i][0];
                try:
                        chromlength = sequences[chrom][1];
                except KeyError:
                        chrom = targetlist[i][0].strip('chr');
                        try: chromlength = sequences[chrom][1];
                        except KeyError: continue;

                GCcontent = 0; bases =0;
                for j in xrange(targetlist[i][1],targetlist[i][2]):
                        bases +=1;
                        if sequences[chrom][0][j] == 'G' or sequences[chrom][0][j] == 'C': GCcontent +=1;
                print targetlist[i][0],targetlist[i][1],targetlist[i][2],
                print targetlist[i][2]-targetlist[i][1],float(GCcontent)/bases,
                SS = sequences[chrom][0][targetlist[i][1]:targetlist[i][2]]; l = len(SS);
                RS = ''.join([revcomp[SS[l-i-1]] for i in xrange(l)]);
                ## print reverse complement of sequence as well 
                #   revsequence = ''.join([revcomp[read[9][l-i-1]] for i in xrange(l)]);
                print SS,RS;


if len(sys.argv)< 3: print 'arguments needed: 1) fastafile 2) region'; sys.exit();
extract_sequences(sys.argv[1],sys.argv[2]);

