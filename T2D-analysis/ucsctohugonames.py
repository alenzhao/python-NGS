#!/usr/bin/python
import os, glob,sys, subprocess,re,random,math, gzip

### written feb 15 2012 
### program that reads a list of genes and a file that maps UCSCids to HUGO/knownIDs and the UCSC transcript file to output transcripts that match the list of genes 
### python ucsctohugonames.py ../129genes ucscgenes-hugonames ucscgenes.ncbi37

### of the 129 genes in current list, all have at least one transcript except C14orf70 

def read_genelist(genelistfile): # list of genes we want to sequence 
        genes = {}; ng =0;
        File = open(genelistfile,'r');
        for line in File: genes[line.strip().split()[0]] = 0;  ng +=1;
        File.close();
        print >>sys.stderr, "read",ng,"genes from file",genelistfile;
#       for gene in genes.iterkeys(): print gene;
        return genes;

# read file that has #kgID displayID and alias (the name in genes{})
def read_genetable(IDfile,genes):
        UCSCids = {}; displayids = {}; # hashtable that store UCSC ids and display ids for genes 
        File = open(IDfile,'r');
        for line in File:
                if line[0] == '#': continue;
                gene = line.strip().split();
                if gene[2] in genes:
                        #genes[gene[2]].append([gene[0],gene[1]]);
                        UCSCids[gene[0].split('.')[0]] = gene[2];
                        displayids[gene[1]] = gene[2];
                        #print gene;
                #### adding _HUMAN was causing genes on different chromosomes to map to the same name, so removed this
                #### 
                """
                elif gene[2].rstrip('_HUMAN') in genes: 
                        #genes[gene[2].rstrip('_HUMAN')].append([gene[0],gene[1]]);
                        UCSCids[gene[0].split('.')[0]] = gene[2].rstrip('_HUMAN');
                        displayids[gene[1]] = gene[2].rstrip('_HUMAN');
                        #print gene;
                """
        return [UCSCids,displayids];

def ucscgenestable(ucscfile,UCSCids,displayids,genes):
        File = open(ucscfile,'r');
        for line in File:
                if line[0] == '#': continue;
                gene = line.strip().split('\t');
                ucscid = gene[0].split('.')[0];  # remove the leading .[123] since that results in few matches
                if ucscid in UCSCids:
                        genes[UCSCids[ucscid]] +=1;
                        print UCSCids[ucscid],line,
                elif gene[10] != '' and gene[10] in displayids:
                        genes[displayids[gene[10]]] +=1;
                        print displayids[gene[10]],line,
                #else:                  print 'missing',line,
                #print gene[0],gene[10];


genes = read_genelist(sys.argv[1]);
[UCSCids,displayids] = read_genetable(sys.argv[2],genes);
ucscgenestable(sys.argv[3],UCSCids,displayids,genes);


for gene in genes.iterkeys():
        if genes[gene] ==0: print >>sys.stderr, 'TRANSCRIPTS',gene,genes[gene];
#for gene in UCSCids.iterkeys(): print gene,UCSCids[gene];




                               
