#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler
import random

# http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.mannwhitneyu.html#scipy.stats.mannwhitneyu 
#last modified aug 21 2014
## calculate p-value for early onset of cases with mutations in GCK genes compared to general background
## python MWU_earlyonset.py Sample_info/AO.list, pvalue = 4.76828523751e-06

from scipy import stats #mannwhitneyu
import scipy 


X = scipy.zeros(18); 

samples = 0; slist = [];
File = open(sys.argv[1]); 
for line in File: 	s = int(line.strip().split()[0]); samples +=1; slist.append(s); 
File.close();
Y = scipy.zeros(samples);
for i in xrange(samples): Y[i] = slist[i]; 


#Y = scipy.zeros(50); 

X[0] = 21; X[1] = 27; X[2] =28; X[3] = 13; X[4] = 39; X[5] = 61; X[6] = 45; X[7] = 33; X[8] = 19; X[9] = 37; X[10] = 30; X[11]= 28; X[12] = 41; X[13] = 24; X[14] = 29; X[15] = 34; X[16] = 30; X[17] = 32; 

print X; print len(Y) 

[U,prob] = stats.mannwhitneyu(X,Y);
print U,prob,math.log(prob,10)


