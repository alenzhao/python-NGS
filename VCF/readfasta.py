#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler
### functions to read fasta file and generate list of reference sequences ###
### useful for samtopileup, hashtable, etc 

def read_fasta(seqfile):
        sequences = {}; reflist = []; lines =0; strlist = [];
        File = open(seqfile,'r');
        if not File: print >>sys.stderr, 'reference sequence file',seqfile,'not found'; sys.exit();

        for line in File:
                if line[0] == '>':
                        if lines > 0:
                                l =0;
                                for k in xrange(len(strlist)): l += len(strlist[k]);
                                sequences[refname][0] += ''.join(strlist); sequences[refname][1] = l;
                                for i in range(len(strlist)): strlist.pop();
                                print >>sys.stderr, 'added seq',refname,(sequences[refname][1]);
                        refname = line.strip().strip('\r').lstrip('>').split()[0];      reflist.append(refname);
                        sequences[refname] = ['',0];    # refsequence and basecalls 
                else: strlist.append(line.strip().strip('\r').upper()); lines +=1;
        l =0;
        for k in xrange(len(strlist)): l += len(strlist[k]);
        sequences[refname][0] += ''.join(strlist); sequences[refname][1] = l;
        for i in range(len(strlist)): strlist.pop();
        print >>sys.stderr, 'added seq',refname,(sequences[refname][1]);
        File.close();
        return [sequences,reflist];

def read_fasta_new(seqfile,fastaindex):
        sequences = [];  File = open(seqfile,'r'); lines =0; strlist = []; fastaseqs =0;
        for line in File:
                if line[0] == '>':
                        refname = line.strip().split()[0].lstrip('>');
                        sequences.append([refname,'',0]); fastaindex[refname] = fastaseqs; fastaseqs +=1;
                        if lines > 0:
                                print >>sys.stderr, 'length of list',len(strlist),sequences[fastaseqs-1][0];
                                sequences[fastaseqs-1][1] = ''.join(strlist);
                                for i in range(len(strlist)): sequences[fastaseqs-1][2] += len(strlist[i]);
                                for i in range(len(strlist)): strlist.pop();
                else: strlist.append(line.strip().upper()); lines +=1;
                        #       print >>sys.stderr, 'length of list',len(strlist);
        sequences[fastaseqs-1][1] = ''.join(strlist);
        for i in range(len(strlist)): sequences[fastaseqs-1][2] += len(strlist[i]);
        for i in range(len(strlist)): strlist.pop();
        File.close();
        for a in range(min(fastaseqs,10)): print >>sys.stderr, sequences[a][0],len(sequences[a][1]),sequences[a][1][0:50];
        return sequences;
