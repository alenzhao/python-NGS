#!/usr/bin/python
import os, glob,sys;

def calculate_rightshift(chromosome,start,end,maxindel):


	## for each base in exon, examine all deletions of length 1-K for VNTR/MS
	for position in xrange(start,end+1):
		for lengthvar in xrange(6,maxindel+1):

			leftpos = position; rightpos = position+lengthvar; rightshift = 0;
			while chromosome[0][leftpos] == chromosome[0][rightpos]:  leftpos +=1; rightpos +=1; rightshift +=1; 
			if rightshift > lengthvar*2: print position,'length',lengthvar,rightshift,chromosome[0][position:rightpos]; 





