#!/usr/bin/python
import os, glob,sys;
### program to convert indels in VCF file to full-length indels that contain the entire haplotype for the reference and indel ### 
### autohr Vikas Bansal, april 3 2012
from readfasta import read_fasta_new


def calculate_rightshift(chrom,position,ref,alt,sequences,fastaindex):

        if len(ref) > len(alt):
                index = fastaindex[chrom];
                lengthvar = len(ref)-len(alt);
                leftpos = position; rightpos = position+lengthvar; rightshift = 0;
                while sequences[index][1][leftpos] == sequences[index][1][rightpos]:  leftpos +=1; rightpos +=1; rightshift +=1;
                newref = sequences[index][1][position-1:position+lengthvar+rightshift+1];
                newalt = sequences[index][1][position-1]+sequences[index][1][position+lengthvar:position+lengthvar+rightshift+1];
                return [newref,newalt,rightshift];

                #print indel[3],indel[4],
                #print '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(chrom,position,indel[2],newref,newalt,indel[5],indel[6],indel[7],indel[8],indel[9])
        elif len(ref) < len(alt):
                index = fastaindex[chrom];
                lengthvar = len(alt)-len(ref);
                if position + 200 < sequences[index][2]: newseq = alt[1:] + sequences[index][1][position:position+200];
                else: newseq = alt[1:] + sequences[index][1][position:sequences[index][2]];
                leftpos = 0; rightpos = lengthvar; rightshift =0; maxl = len(newseq);
                while newseq[leftpos] == newseq[rightpos] and rightpos < maxl: leftpos +=1; rightpos +=1; rightshift +=1;

                newref = sequences[index][1][position-1:position+rightshift+1];
                newalt = alt + sequences[index][1][position:position+rightshift+1];
                #print 'INS',rightshift,newref,newalt;
                #print line

                #print chrom,position,type,bases, sequences[index][1][position-1-rightshift:position-1],bases,rightshift;
        #       print line, chrom,position, rightshift;
                return [newref,newalt,rightshift];
        else:  # it is a SNP 
                return [ref,alt,0];

def read_indels(indelfile,seqfile):
        fastaindex = {}; sequences = read_fasta_new(seqfile,fastaindex);

        if indelfile == '-' or indelfile == 'stdin' or indelfile == 'sys.stdin': File = sys.stdin;
        else: File = open(indelfile,'r');
        for line in File:
                if line[0] == '#': print line.strip(); continue;

                indel = line.strip('\n').split('\t'); chrom = indel[0]; position = int(indel[1]);

                ## print line for SNPs
                if len(indel[3]) == len(indel[4]): print line.strip(); continue;
                #if len(altalleles) ==1: continue;
                ref = indel[3]; altalleles = indel[4].split(',');

                for i in xrange(len(altalleles)):
                        s1=len(ref)-1; s2 = len(altalleles[i])-1;
                        if s1 ==s2: continue;
                        if s1> s2 and s2 > 0: # both alleles are greater than 1 in length:
                                while s2 > 0 and ref[s1] == altalleles[i][s2]: s1 -=1; s2 -=1;
                                ref0 = ref[0:s1+1]; alt0 = altalleles[i][0:s2+1];
                        elif s2> s1 and s1 > 0: # both alleles are greater than 1 in length:
                                while s1 > 0 and ref[s1] == altalleles[i][s2]: s1 -=1; s2 -=1;
                                ref0 = ref[0:s1+1]; alt0 = altalleles[i][0:s2+1];
                        else: ref0 = ref; alt0 = altalleles[i];

                        [newref,newalt,shift] = calculate_rightshift(chrom,position,ref0,alt0,sequences,fastaindex);
                        #print indel[3],indel[4],ref0,alt0,shift,'|',
                        print '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(chrom,position,indel[2],newref,newalt,indel[5],indel[6],indel[7],indel[8],indel[9])
                        #continue;
                ## take care of case where there are multiple alternate alleles...
        if indelfile == '-' or indelfile == 'stdin' or indelfile == 'sys.stdin': pass;
        else: File.close();
#if len(sys.argv) < 3: print  'python leftjustify-indels.py indel-file reference.fasta FORMAT'; sys.exit();
#if len(sys.argv) > 3: FORMAT = sys.argv[3];
#read_indels(sys.argv[1],sys.argv[2]);
