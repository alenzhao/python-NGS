#! /usr/bin/env python
import sys, os, glob, string

#### compare two VCF files and calculate # of variants that differ ####
#### printflag = 1, prints the variants that are not in File 2 but in File1 1
def compareVCF(vfile1,vfile2,printflag):
        File1 = open(vfile1,'r');
        variants = {};
        for line in File1:
                if line[0] == '#': continue;
                var = line.strip().split();
                variants[var[0]+':'+var[1]+':'+var[3]+':'+var[4]] = [line.strip(),1];
        File1.close();

        shared = 0; unique1 = 0; unique2=0;
        File2 = open(vfile2,'r');
        for line in File2:
                if line[0] == '#': continue;
                var = line.strip().split('\t');
                try:
                        match = variants[var[0]+':'+var[1]+':'+var[3]+':'+var[4]];
                        match[1] = 0; # shared 
                except KeyError:
                        variants[var[0]+':'+var[1]+':'+var[3]+':'+var[4]] = [line.strip(),2];
        File2.close();

        nonshared = [];
        for var in variants.iterkeys():
                if variants[var][1] == 1:
                        unique1 +=1;
                        if printflag ==1: print 'FILE1',variants[var][0];
                elif variants[var][1] == 2:
                        unique2 +=1;
                        if printflag ==2: print 'FILE2',variants[var][0];
                elif variants[var][1] == 0: shared +=1;
        print 'variants shared:',shared,'unique to VCF1:',unique1,'unique to VCF2:',unique2;
~                                                                                                                                                                                              
~                                                                     
