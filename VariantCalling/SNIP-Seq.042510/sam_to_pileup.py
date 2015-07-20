#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler
from array import *

######################################################################

bti = {'A':0,'C':1,'G':2,'T':3,'a':0,'c':1,'g':2,'t':3,'N':4};
itb = {0:'A', 1:'C',2:'G',3:'T'};

def flagtostring(flag,fstring):
	for i in range(10):  fstring[i] = `(flag%2)`; flag /= 2; 

#a = array('c','000000000000'); flagtostring(147,a); print a;  sys.exit();


def read_fasta(seqfile):
	sequences = {}; reflist = []; lines =0; strlist = []; 
	File = open(seqfile,'r');
	if not File: print >>sys.stderr, 'reference sequence file',seqfile,'not found'; sys.exit();

	for line in File:
		if '>' in line: 
			if lines > 0:
				sequences[refname][0] += ''.join(strlist); 
				print >>sys.stderr, 'added seq',refname,len(sequences[refname][0]);
				for i in range(len(strlist)): strlist.pop();
			refname = line.strip().strip('\r').lstrip('>').split()[0];
			sequences[refname] = ['',{},{}];    # refsequence and basecalls 
			reflist.append(refname);
		else: strlist.append(line.strip().strip('\r').upper()); lines +=1; 
	sequences[refname][0] += ''.join(strlist); 
	print >>sys.stderr, 'added seq',refname,len(sequences[refname][0]);
	File.close();
	return [sequences,reflist];
	

def print_pileup(reflist,current,index,start,end):
	for j in range(start,end):
		bases = '@'; qvalues = '@'; mvalues = '@'; positions = ''; mismatches = ''; coverage =0; baseslist = [];
		if j in index[1]: 
			T = index[1][j]; # T is the string that has basecalls for position j in the refsequnece 
			coverage = (len(T)-1)/6;
			qvalues += ''.join([T[6*k+2] for k in xrange(coverage)]);
			positions  += ''.join([str(ord(T[6*k+3])-33) + ',' for k in xrange(coverage)]);
			mvalues += ''.join([T[6*k+4] for k in xrange(coverage)]);
			mismatches += ''.join([T[6*k+6] for k in xrange(coverage)]);
			for k in range(coverage): 
				strand  = T[6*k+5];
				if strand == '+' and index[0][j] == T[6*k+1]: baseslist.append(',');
				elif strand == '-' and index[0][j] == T[6*k+1]: baseslist.append('.');
				elif strand == '+': baseslist.append(T[6*k+1]);
				else: baseslist.append(T[6*k+1].lower()); 
			bases += ''.join([baseslist[k] for k in xrange(coverage)]);

			print reflist[current],j+1,index[0][j],coverage,bases,qvalues,mvalues,positions,mismatches,
			if j in index[2]: print ''.join(index[2][j]); del index[2][j];
			else: print;
			del index[1][j]; 
		else:
			print reflist[current],j+1,index[0][j],coverage,bases,qvalues,mvalues,positions,mismatches; 
			pass;


def samtopileup(readfile,seqfile):
	[sequences,reflist] = read_fasta(seqfile); 
	#for a,b in sequences.iteritems(): refsequence = b[0]; 	print >>sys.stderr, a,len(b[0]);
	current =0 ; index = sequences[reflist[current]]; last =0;

	File=open(readfile,'r');  reads =0;
	if not File: print >>sys.stderr, 'samfile',readfile,'not found'; sys.exit();
	editstring = array('c'); 
	for i in range(1000): editstring.append('.');

	for s in File:
		if s[0] == '@': continue;
		line = s.strip().split();
		if line[2] == '*': continue;
		reads +=1; 
		if reads%100000 ==0: print >>sys.stderr, 'processed',reads;
		if len(line) < 11: continue;
		readid = line[0]; flag = int(line[1]); locus = line[2]; start = int(line[3]); mq = int(line[4]); cigarstring = line[5]; 
		read = line[9]; quality = line[10]; 
		try: refsequence = sequences[locus];  # index into refnames table 
		except KeyError: print >>sys.stderr, 'could not find sequence name',locus,'in fasta file'; continue; 
		offset =0;

		if flag & 4 == 4: continue;  # unmapped read 
		if flag & 16 == 16: strand = '-';
		else: strand = '+';
		if flag & 1 == 1:
			if flag & 64 == 64: read12 = 1;
			elif flag & 128 == 128 : read12 = 2;
			else: read12 = 0;
		else: read12 = 0;
		
		if start -offset -1 < 0 or start -offset -1 + len(read) >= len(refsequence[0]): continue; 
		if cigarstring == "*" or 'H' in cigarstring or 'P' in cigarstring or 'N' in cigarstring: continue;

		readlength = len(read);
	
		l1 = 0; e = 0; #editstring = '';
		for k in xrange(len(cigarstring)):
			if cigarstring[k] == 'M': 
				for f in xrange(int(cigarstring[l1:k])): editstring[e] = '.'; e +=1;
				l1 = k+1;
			if cigarstring[k] == 'S': 
				for f in xrange(int(cigarstring[l1:k])): editstring[e] = 's'; e +=1;
				l1 = k+1;
			elif cigarstring[k] == 'D': 
				for f in xrange(int(cigarstring[l1:k])): editstring[e] = '-'; e +=1;
				l1 = k+1;
			elif cigarstring[k] == 'I': 
				for f in xrange(int(cigarstring[l1:k])): editstring[e] = '+'; e +=1;
				l1 = k+1;

		if read12 ==2: delta = readlength;
		else: delta =0;

		l1 = 0; l2 = 0;		mismatches =0;
		for i in xrange(e):
			if editstring[i] == '.':
				if read[l1].upper() != refsequence[0][start+l2-offset-1]: mismatches +=1;
				l1 +=1; l2 +=1;
			elif editstring[i] == 's':
				l1 +=1; l2 +=1;
			elif editstring[i] == '-':
				if editstring[i+1] != '-' : mismatches +=1;
				l2 +=1;
			elif editstring[i] == '+':
				if editstring[i+1] != '+' : mismatches +=1;
				l1 += 1;  

	#	if 'I' in cigarstring or 'D' in cigarstring or 'M' in cigarstring: 
	#		print read,start-offset-1,editstring,cigarstring,mismatches;
	#		print refsequence[0][start+offset-1:start+offset-1+readlength];
		
		l1 = 0; l2 = 0; flag = 0; bases = '';
		for i in xrange(e):
			if editstring[i] == '.':
				flag =0;
				if strand == '+': psb = chr(l1+33+1+delta);
				else: psb = chr(len(read) -l1 + 33+delta);
				mqb = chr(min(123,mq+33)); 
				try: refsequence[1][start+l2-offset-1].fromstring(read[l1]+quality[l1]+ psb + mqb + strand + chr(mismatches+48));
				except KeyError: refsequence[1][start+l2-offset-1] = array('c','_' + read[l1]+quality[l1]+ psb + mqb + strand + chr(mismatches+48));
				l1 +=1; l2 +=1;
			elif editstring[i] == 's':
				l1 +=1; l2 +=1;
			elif editstring[i] == '-':
				if flag ==0: bases = refsequence[0][start+l2-offset-1]; flag =1;
				else: bases += refsequence[0][start+l2-offset-1]; 
				if editstring[i+1] != '-': 
					if strand == '+': 
						try: refsequence[2][start+l2-offset-1].fromstring('-' + bases);
						except KeyError: refsequence[2][start+l2-offset-1] = array('c','-' + bases);
					else:
						try: refsequence[2][start+l2-offset-1].fromstring('-' + bases.lower());
						except KeyError: refsequence[2][start+l2-offset-1] = array('c','-' + bases.lower());
				l2 +=1;
			elif editstring[i] == '+':
				if flag ==0: bases = read[l1]; flag =1;
				else: bases += read[l1];
				if editstring[i+1] != '+' : 
					if strand == '+': 
						try: refsequence[2][start+l2-offset-1].fromstring('+' + bases);
						except KeyError: refsequence[2][start+l2-offset-1] = array('c','+' + bases);
					else:
						try: refsequence[2][start+l2-offset-1].fromstring('+' + bases.lower());
						except KeyError: refsequence[2][start+l2-offset-1] = array('c','+' + bases.lower());
				l1 += 1;  
		
		while locus != reflist[current]: 
			# print pileup from previous to end of len(index[0])
			print_pileup(reflist,current,index,last,len(index[0])); last = 0;
			current +=1; 
			if current < len(reflist): index = sequences[reflist[current]];
			else: break; 
		if current >= len(reflist): break; 
		print_pileup(reflist,current,index,last,start-1); last = start-1;

	print_pileup(reflist,current,index,last,len(index[0])); 
	last =0;
	current +=1; 
	while current < len(reflist):
		index = sequences[reflist[current]];
		print_pileup(reflist,current,index,last,len(index[0])); 
		current +=1; 
		


	File.close();


if len(sys.argv)< 3: print 'python sam_pileup.py SAMfile referencesequence.fasta'; sys.exit();
samtopileup(sys.argv[1],sys.argv[2]);


