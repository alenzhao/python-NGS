#!/usr/bin/python
import os, glob,sys, subprocess,re,random,math, gzip

## code to extract exon start-end coordinates for genes in a list using ucsc annotations/refseq
## BUG FIXED in overlap function on march 21 2013, caused KCNJ11 exon piece to be missed, check if other exons missed 

def compute_exonlist(ucscfile):
        GENES = {};
        File= open(ucscfile,'r');
        for line in File:
                if line[0] == '#': continue; ## header line
                transcript = line.strip().split(); genename = transcript[0];
                if genename not in GENES: GENES[genename] = [0];

                exonstarts = transcript[9].rstrip(',').split(','); exonends = transcript[10].rstrip(',').split(',');
                if GENES[genename][0] ==0:
                        GENES[genename].append([]);
                        GENES[genename].append([]);
                        GENES[genename].append([]);
                        GENES[genename].append([]);
                        GENES[genename].append(int(transcript[8]));

                GENES[genename][0] +=1;
                #print genename,GENES[genename][0],exonstarts,exonends;

                ## add the exon to the list of exons if it is not already present
                for s in xrange(len(exonstarts)):
                        if (int(exonstarts[s]),int(exonends[s])) not in GENES[genename][3]: GENES[genename][3].append((int(exonstarts[s]),int(exonends[s])));

                # add (txstart,txend) pair to list of transcripts 
                if (int(transcript[4]),int(transcript[5])) not in GENES[genename][4]:
                        GENES[genename][4].append((int(transcript[4]),int(transcript[5])));
                if int(transcript[8]) > GENES[genename][5]: GENES[genename][5] = int(transcript[8]);

                GENES[genename][1].append(transcript[2]);
                ## add strand information 
                GENES[genename][2].append(transcript[3]);
                #print genename,transcript;

        TOTALBASES =0; TOTALSPAN =0; PROMBASES =0; PRINTEXON =1; FLANKING = 1000;

        # need to do overlap-merge on sorted list of exons for each gene...
        for a,b in GENES.iteritems():
                if len(b) < 3: print a,'MISSING';
                else:
                        b[1].sort();
                        for i in xrange(len(b[1])-1):
                                if b[1][i] != b[1][i+1]: print 'different chromosomes';
                        b[3].sort();

                        #print b[3];
                        totalexonlength = 0; span =0; prev = b[3][0][0];
                        for i in xrange(len(b[3])-1):
                                if b[3][i][1] < b[3][i+1][0]:
                                        totalexonlength += b[3][i][1]-prev+1;
                                        if PRINTEXON: print b[1][0],prev,b[3][i][1],a,b[2][0],'EXON';
                                        prev = b[3][i+1][0];
                                elif b[3][i][1] > b[3][i+1][1]: # this exon encompasses the next exon completely 
                                        print 'merging next exon',a,b[3][i][0],b[3][i][1],b[3][i+1][0],b[3][i+1][1];
                                        b[3][i+1] = b[3][i]; #prev = b[3][i+1][0]; 
                                        pass; # bug fixed here march 21 2013
                                        #if PRINTEXON: print b[1][0],prev,b[3][i][1],a,b[2][0],'EXON_OVERLAPPING';
                                        #prev = b[3][i+1][0];   
                                else:
                                        pass;

                        totalexonlength += b[3][len(b[3])-1][1]-prev+1;
                        if PRINTEXON: print b[1][0],prev,b[3][len(b[3])-1][1],a,b[2][0],'EXON';

                        span = b[4][0][1]-b[4][0][0]+1;
                        #for s in b[3]:         print b[1],s[0],s[1],a; 
                        TOTALBASES += totalexonlength;
                        TOTALSPAN += span;

                        for i in xrange(len(b[4])):
                                if b[2][0] == '+': print b[1][0],b[4][i][0]-FLANKING,b[4][i][0],a,b[2][0],'PROMOTER';
                                elif b[2][0] == '-': print b[1][0],b[4][i][1],b[4][i][1]+FLANKING,a,b[2][0],'PROMOTER';
                                PROMBASES += FLANKING;
                                #TOTALBASES += FLANKING;

                        print 'GENESUMMARY',a,b[0],b[1][0],b[5],len(b[3]),totalexonlength,'SPAN',span,b[2][0],b[4];

        print >>sys.stderr, 'TOTALBASES',TOTALBASES,PROMBASES,TOTALSPAN;

        File.close();


#####################################################################################################

compute_exonlist(sys.argv[1]);


