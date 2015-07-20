#! /usr/bin/env python
import sys, os, glob, string, math, array
import re
space = re.compile(r'\s+');
import subprocess;

#cmd = 'ls -l '; os.system(cmd); sys.exit();


from shutil import copy;


def int_to_kmer(kmerint,kmerlength):
	s ='';	
	for i in range(kmerlength):
		s += int_to_base(kmerint %4);
		kmerint = kmerint /4;
	slist = list(s); slist.reverse(); slist = ''.join(slist); 
	return slist;
	

def int_to_base(b):
        if b == 0: return 'A';
        if b == 1: return 'C';
        if b == 2: return 'G';
        if b == 3: return 'T';
        else: return '-';


def base_to_int(b):
        if b == 'A' or b == 'a': return 0;
        elif b == 'C' or b == 'c': return 1;
        elif b == 'G' or b == 'g': return 2;
        elif b == 'T' or b == 't': return 3;
        elif b == 'N' or b == 't': return 0;
	else: return 0; 

kmerlength = 10;
if len(sys.argv) >2: kmerlength  = int(sys.argv[2]);


#a = 'vikas'; a1 = list(a); a1.reverse(); a1 = ''.join(a1); print a1; sys.exit();

subsvec = [0,0,0,0]; 
subsvec[0] = base_to_int('A'); subsvec[1] = base_to_int('C'); subsvec[2] = base_to_int('G'); subsvec[3] = base_to_int('T');
for i in range(kmerlength-1): subsvec[0] = subsvec[0] << 2; subsvec[1] = subsvec[1] << 2; subsvec[2] = subsvec[2] << 2; subsvec[3] = subsvec[3] << 2;
print subsvec[0], subsvec[1], subsvec[2], subsvec[3];




#subprocess.call("sleep 30", shell= True );



if len(sys.argv)> 1: seqfile = sys.argv[1];
else: print " enter filename to read ...."; sys.exit();
if seqfile.find('/') >=0: dirname = seqfile; directory =1;
else: directory = 0; 


hashtablesize = 1; 
for i in range(kmerlength): hashtablesize *= 4;
kmercounts = [];
for i in range(hashtablesize):	kmercounts.append(0);




listofgenomefiles = [];
if directory ==1: 
	for file in glob.glob(dirname + "/*.fa"): listofgenomefiles.append(file);
else: listofgenomefiles.append(seqfile);

print listofgenomefiles;


for seqfile in listofgenomefiles:

	kmerstring = [];
	kmerint = 0;
	File = open(seqfile,'r'); 
	lines=0;
	while 1:
		line = File.readline(); lines +=1;
		if not line: break;
	File.close();
	totalines = lines;

	stdstring ='read file ' + seqfile + ' with ' + str(totalines) + ' lines \n'; 
#	sys.stderr.write(stdstring);
 
	File = open(seqfile,'r'); 
	lines =0; wrap =0;kmers =0;
	while 1:
		line = File.readline().rstrip('\n'); lines +=1;
		if not line: break;
		if line.find('>') >= 0: pass; 
		else:  
			#print line.rstrip('\n'); 
			for i in range(len(line)-1): 
				kmerint = kmerint << 2; kmerint += base_to_int(line[i]); 
#				if line[i] == 'A' or line[i] == 'a' or line[i] == 'c' or line[i] =='C' or line[i] =='g' or line[i] == 'G' or line[i] == 't' or line[i] == 'T':
					#if wrap <= 0: kmercounts[kmerint] += 1;
					#else: wrap -= 1;
				#else: wrap = kmerlength;

				kmercounts[kmerint] += 1; kmers +=1;
				kmerstring.append(line[i]);
				#print kmerstring,kmerint;
				if len(kmerstring) == kmerlength : 
			 		kmerint -= subsvec[base_to_int(kmerstring[0])];
				 	#print subsvec[base_to_int(kmerstring[0])];
					kmerstring.pop(0); 

		if lines % int(totalines/100) ==0: stdstring = 'done for ' + str(int(float(lines*100)/totalines)) + ' percent of file\n'; 
		#sys.stderr.write(stdstring); 
	#	if lines > 10000: break;

	File.close();

	for i in range(len(kmerstring)): kmerstring.pop(0);
	#sys.stderr.write('kmer count '); sys.stderr.write(str(kmers)); sys.stderr.write('\n'); 

	outputfile = seqfile + '.' + str(kmerlength)+ '.hashtable'; 
	Fileoutput = open(outputfile,'w');
	print 'outputting hash table kmer count.... ',seqfile,outputfile;
	Fileoutput.write('outputting hash table kmer count.... \n' + seqfile);
	for i in range (hashtablesize): Fileoutput.write(str(i) + "\t" + str(int_to_kmer(i,kmerlength)) + "\t" + str(kmercounts[i]) + '\n'); 
	Fileoutput.close();
	subprocess.call("sort -g -k 3 " + outputfile + " > " + outputfile + ".sorted", shell=True); 


#python count-kmers.py chr2.fa 10 > chr2.fa.kmercounts_k10N ;  more chr2.fa.kmercounts_k10N | sort -g -k 3  | awk 'BEGIN { c=0; k =0; } { if ($3 == c) k++;  else {print c,k; k=1;  } c =$3;  } '  > counts10_chr2

print 'outputting hash table kmer count.... \n';
totalkmercounts =0; 
for i in range (hashtablesize): print i,int_to_kmer(i,kmerlength),kmercounts[i]; totalkmercounts += kmercounts[i];
print "total kmers hashed: ",totalkmercounts; 

sys.exit();







