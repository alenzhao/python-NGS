#!/usr/bin/python
import os, glob,sys, subprocess,re,random,math, gzip

## post-process GSNAP SAM/BAM file to remove reads with multiple alignments on independent lines (read-name sorted)

def find_bestpair(alignments,a):
        if a ==2: return [0,1];
        if a ==3: return 1;


def process_gsnap(gsnapfile):
        if gsnapfile == "stdin" or gsnapfile == "-": File = sys.stdin;
        else: File = open(gsnapfile,'r');
        alignments = []; readid = ""; a =0; a1 = 0; a2 = 0;
        for line in File:
                read = line.split(); flag = int(read[1]) & 64;  mq = int(read[4]);
                if read[0] == readid or readid == "": alignments.append([flag,mq,line]); readid = read[0]; a +=1;
                else:
                        a1 = 0; a2= 0;
                        for i in xrange(a):
                                if alignments[i][0] == 64: a1 +=1;
                                else: a2 +=1;
                        if a == 2 and (alignments[0][1] >= 20 or alignments[1][1]>= 20) :  print alignments[0][2],alignments[1][2],
                        """
                        elif a1 > 1 and a2 > 1: pass: 
                        if a > 2: 
                                for i in xrange(a): print alignments[i][0],alignments[i][1],alignments[i][2], 
                                print a1,a2;
                        """

                        for i in xrange(a): alignments.pop();
                        a = 1; mq = int(read[4]);
                        alignments.append([flag,mq,line]); readid = read[0]; a = 1;





process_gsnap(sys.argv[1]);
