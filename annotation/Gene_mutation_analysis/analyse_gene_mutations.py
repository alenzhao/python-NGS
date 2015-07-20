import sys,os, math
from optparse import OptionParser

import multinomial_test
import read_ExAc

## within each exon, mutations can further be partitionned based on conservation/PolyPhen score and analysis done with higher resolution..

## estimate gene mutation rate per bp using silent mutations -> estimate expected # of missense mutations...

## if mutation_type == 'S/NS', program compares mutation counts between silent/missense per exon 

revcomp = {'A':'T', 'T':'A', 'C':'G', 'G':'C','a':'T', 't':'A', 'c':'G', 'g':'C', 'N':'N','n':'N' }
VERBOSE =0;

#########################################################################################

### python xx.py file_all_mutations disease_mutations/VCF from ExAc type_of_file type_of_mutations

## open file with set of all stop/missense mutations in a single gene | format = Exon 2 70 24 GAG TAG G T E * trimer TGA T 20 127 1.80084386596e-09 position 
if len(sys.argv) > 5 and sys.argv[5] == "verbose": VERBOSE =1

if len(sys.argv) > 4: mutation_type = sys.argv[4]; print >>sys.stderr, "will only analyze mutations of type", mutation_type;
else: mutation_type = "None"; 


#### read file with list of all possible mutations in the coding sequence of gene 
###########################################################################################################
allmutations= sys.argv[1]; 
mutlist = {}; exonlist = {}; mutlist_position = {}; STRAND = '+'; 
first = -1; last = -1; nmuts =0; chromosome = "None";
File = open(allmutations,'r'); 
for line in File:
	mut = line.strip().split();
	if mutation_type != 'S/NS' and mutation_type != 'None' and mut[18] != mutation_type: continue; 
	if mutation_type == 'S/NS' and mut[18] != 'silent' and mut[18] != 'missense': 
		continue; ## compare 

	exon = int(mut[1]); AApos = int(mut[3]); mutprob = float(mut[15]); nmuts +=1;
	
	## variant info chr:position:ref:alt
	variant = mut[16].split(':'); position = int(variant[1]); 
	if variant[2] != mut[6] and STRAND == '+': STRAND = '-';
	if chromosome == "None": chromosome = variant[0]; 
	if position < first or first == -1: first = position; 
	if position > last: last = position; 


	mutlist_position[(variant[0],position,variant[2],variant[3])] = [exon,mutprob,mut[17],mut[18]]; ## useful for VCF data comparison 

	## we can have two missense mutations at same position, this is only used to find exon corresponding to each mutation  
	try: 	mutlist[AApos][2] += 1; 
	except KeyError: mutlist[AApos] = [exon,mut[8],1]; 

	try: 	exonlist[exon][0] +=1; 
	except KeyError: exonlist[exon] = [1,0.0,0,0,0.0];

	if mutation_type != 'S/NS': exonlist[exon][1] += mutprob;
	else:
		if mut[18] == 'silent': exonlist[exon][1] += mutprob;
		elif mut[18] == 'missense': exonlist[exon][4] += mutprob;

File.close(); 
print >>sys.stderr, "read",nmuts,"mutations from file",allmutations,'firstpos:',first,'lastpos:',last; 

###########################################################################################################

if sys.argv[3] == 'controls': 
	## open ExAC VCF file or any other VCF file and update counts of observed mutations per exon  
	mutations = {}; sharedmuts = 0;
	variant_list = read_ExAc.read_vcf(sys.argv[2],chromosome,first,last,mutations);
	for variant in variant_list: 
		try: 
			annotation = mutlist_position[(variant[0],variant[1],variant[2],variant[3])]; 
			exon = annotation[0]; 
			if mutation_type != 'S/NS': exonlist[exon][2] +=1; 
			else: 
				if annotation[3] == 'silent': exonlist[exon][2] +=1; 
				elif annotation[3] == 'missense': exonlist[exon][3] +=1; 
				
			sharedmuts+=1;
			if VERBOSE ==1: print >>sys.stderr, variant,annotation,exon; 
		except KeyError: pass; 

	print >>sys.stderr, "shared_mutations in VCF file ", sharedmuts;
#sys.exit();

###########################################################################################################

elif sys.argv[3] == 'disease': 

	## open file with set of observed disease stop mutations 
	## p.Y1944X c.5832C>A 27 Bamshad2011 format of disease mutation file
	diseasemutations = sys.argv[2]; 
	File = open(diseasemutations,'r');
	for line in File:
		mut = line.strip().split();
		if 'missense' in mut[0]: 
			AApos = int(mut[2])
		elif 'X' in mut[0]: 
			AAchange = mut[0].lstrip('p.').rstrip('X').rstrip('*'); AA = AAchange[0]; AApos = int(AAchange[1:]);
		elif 'M' in mut[0]: 
			AAchange = mut[0].lstrip('p.').rstrip('M').rstrip('*'); AA = AAchange[0]; AApos = int(AAchange[1:]);
		else: continue;
		#print mut,AApos;
		try: 
			exon = mutlist[AApos]; 
			#print 'mutation found',exon,mut;
			exonlist[exon[0]][2] +=1; 
		except KeyError: 
			print 'mutation not found in list',mut; 
	File.close();

## file with mutations, chr7:4444:- is colum 6 |  NM_000162.3:c.445A-C is column 7 
elif sys.argv[3] == 'database':
	diseasemutations = sys.argv[2]; sharedmuts = 0;
	File = open(diseasemutations,'r');
	for line in File:
		mut = line.split('\t');
		if len(mut) < 12 or mut[5] == 'null': continue;
		alleles = mut[11].split('[')[1].split(']')[0].split('/');
		chromosome = mut[5].split(':')[0].strip('chr'); position = int(mut[5].split(':')[1])
		if STRAND == '+': ref = alleles[0]; alt = alleles[1]; 
		if STRAND == '-': 
			ref = ''.join([revcomp.get(alleles[0][i],'X') for i in range(len(alleles[0]))])
			alt = ''.join([revcomp.get(alleles[1][i],'X') for i in range(len(alleles[1]))])

		print mut[5],mut[6],mut[7],chromosome + ':' + `position` + ':' + ref + ':' + alt;
		try: 
		 	annotation =  mutlist_position[(chromosome,position,ref,alt)] 
			exon = annotation[0]; exonlist[exon][2] +=1; 
			sharedmuts+=1;
			#print >>sys.stderr, chromosome,position,ref,alt,annotation; 

		except KeyError: pass;




#mutrate= 8e-8;
pvals = []; counts = []; mutations = 0; pvalsum =0;
for e in exonlist.iterkeys(): 
	if exonlist[e][1] > 1e-100: 
		if exonlist[e][2] >=0: print 'EXON',e,exonlist[e][0],exonlist[e][1],exonlist[e][2],exonlist[e][4],exonlist[e][3],
		#print 'mutation-rate:',exonlist[e][1]/(exonlist[e][2]+0.0001),exonlist[e][1]/mutrate;
		if mutation_type == 'S/NS': # compare rate of missense to silent mutations per exon 
			R = exonlist[e][2] + exonlist[e][3]; 
			A = exonlist[e][2];
			ep = exonlist[e][1]/(exonlist[e][1] + exonlist[e][4]); 
			[pl,pr] = multinomial_test.binomial_pvalue(R,A,ep); 
			OR = exonlist[e][2]*exonlist[e][4]/(exonlist[e][1]*exonlist[e][3]+0.000001);
			print pl,pr,OR;
		else: print;

		pvals.append(exonlist[e][1]); counts.append(exonlist[e][2]);
		mutations += exonlist[e][2]; pvalsum += exonlist[e][1];

if len(pvals) > 1: ## at least two exons
	for p in xrange(len(pvals)): 
		pvals[p] /= pvalsum; 
		pvalue = [0,0];
		pvalue = multinomial_test.binomial_pvalue(mutations,counts[p],pvals[p]);
		if pvalue[0] < -2: print 'Exon_sig-',p+1,p,pvals[p],counts[p],'expected',int(pvals[p]*mutations),pvalue[0],pvalue[1];
		elif pvalue[1] < -2: print 'Exon_sig+',p+1,p,pvals[p],counts[p],'expected',int(pvals[p]*mutations),pvalue[0],pvalue[1];
		print 'Exon',p+1,p,pvals[p],counts[p],'expected',int(pvals[p]*mutations),pvalue[0],pvalue[1];

#multinomial_test.single_exon_test(pvals,counts);
#multinomial_test.multinomial_test(pvals,counts);

#counts[12] += 10; mutations += 10;




