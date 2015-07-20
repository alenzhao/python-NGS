#!/usr/bin/python
import os, glob,sys, subprocess,re,math,random 

epsilon = 1e-8;
## first implemented jan 27 2015, vikas bansal, vibansal@ucsd.edu 
## KCNB1 G379R, score = 1.0 S347R, score = 1.0  T374I , score =1.0

"""
1. read gene file from dbSNF and store list of all possible non-syn mutations in the gene... what about synonymous mutations ??
2. also get first and last position of gene and use tabix to get list of observed mutations in population data VCF+tabix 
3. read mutation matrix file and store as a table.... fixed 
4. 

how to get information about all possible synonymous mutations in the gene: ?? need transcript model...

 python pvalue.py /media/drive2/Variant-call-datasets/dbNSFPv2.8.baylor/KCNB1.data ~/Public/tools/reference-genomes/NCBI37/human_g1k_v37.fasta

"""

from read_fasta_bychrom import make_fasta_index, read_chromosome

def generate_syn_mutations(CHROM,position,codon,strand): 
	## if strand is negative, codon is backward... TGA is actually TCA on + strand 
	## what about codon that spans two exons.... position will be tricky....
	for i in xrange(3): print 1;



mut_table = {}; bases = ['A','C','G','T']
File = open("trimer_mu_matrix.lis",'r');
for line in File:
	if line[0] == '#': continue;
	trimer = line.split();
	for i in xrange(4): mut_table[(trimer[0],bases[i])] = float(trimer[i+1]); 
File.close();
#for trimer in mut_table.keys(): print trimer,mut_table[trimer]




## fasta_file 
offset_index = make_fasta_index(sys.argv[2]); current_chrom = '-';

if '/' in sys.argv[1]: dbSNFP_file = sys.argv[1]; 
else:
	gene_name = sys.argv[1]; 
	dbSNFP_file = 'GENE_DATA/' + gene_name + '.data';
	if not os.path.isfile(dbSNFP_file):
		# find chromosome for the gene 
		refseqfile = open('/home/vbansal/CODE/JOINTCODE-coral/PYTHON-scripts/VCF-programs/refseq.genes.ncbi37','r');
		for line in refseqfile: 
			tx = line.split(); 
			#print tx[12]
			if tx[12] == gene_name: chromosome = tx[2]; print chromosome, gene_name; break; 

		## read the chromosome file for that gene 
		dirpath = "/media/drive2/Variant-call-datasets/dbNSFPv2.8.baylor/dbNSFP2.8_variant." + chromosome; 
		f = open(dbSNFP_file,'wb');
		print dirpath,gene_name 
		#subprocess.call(tabix_command,stdout=f);
		subprocess.call(["grep",gene_name,dirpath],stdout=f);
		f.close();
	#sys.exit()


################################ read gene file from dbNSFP ######################################
mutations = {}; TXS = {}; firstpos= -1; lastpos = -1; mlist = {}; prevpos = 0; currentpos = -1;
syns = 0; missense = 0; stop = 0; c = 0; 
File = open(dbSNFP_file,'r'); 
for line in File: 
	var = line.split('\t');
	print var[0],var[1],var[2],var[3],var[4],var[5],var[10],var[22],var[23],var[24],var[25],var[29],var[30],var[31]
	if line[0] == '#': continue;
	chrom = var[0]; pos = int(var[1])-1; currentpos = int(var[1]);
	strand = var[15]; codon = var[16]; 
	try: 
		offset = int(var[18]); 
		#if strand == '+' and offset ==3 and currentpos != prevpos: print strand,codon,offset,currentpos-2; 
		#elif strand == '-' and offset ==1 and currentpos != prevpos: print strand,codon,offset,currentpos-2; 
	except ValueError: pass;
	#print strand,codon,offset;

	if pos < firstpos or firstpos < 0: firstpos = pos; 
	if pos > lastpos or lastpos < 0: lastpos = pos+2; 

	if chrom != current_chrom:
		if chrom in offset_index: CHROM  = read_chromosome(sys.argv[2],offset_index,chrom); 
		else: continue;
		current_chrom = chrom;

	if len(var) < 31: continue;
	transcripts= var[22].split(';'); positions = var[23].split(';'); 
	if len(positions) != len(transcripts) and len(positions) > 0: 
		print positions,transcripts,line
		continue;
	mut_rate = mut_table[(CHROM[0][pos-1:pos+2],var[3])];
	if var[5] == 'X': print 'stopgain',var[0],var[1],var[2],var[3],var[4],var[5],var[10],var[22],var[23],var[24],var[25],var[29],var[30],var[31]
	mutations[(var[0],int(var[1]),var[2],var[3])] = [var[4],var[5],var[10],var[29],var[30],CHROM[0][pos-1:pos+2],[],mut_rate,0,[],[]]; 

	for t in xrange(len(transcripts)):
		#print var[0],pos,CHROM[0][pos-2:pos+3],var[2],var[3],var[4],var[5]
		mutations[(var[0],int(var[1]),var[2],var[3])][6].append(transcripts[t]);  
		mutations[(var[0],int(var[1]),var[2],var[3])][9].append(int(positions[t]));  
		#mutations[(var[0],int(var[1]),var[2],var[3])][6].append([transcripts[t],int(positions[t])]);  
		try: TXS[transcripts[t]] += 1; 
		except KeyError: TXS[transcripts[t]] = 1; 
		if var[5] == 'X': stop +=1; 
		else: missense +=1; 

	prevpos = currentpos; 
	
File.close();

print 'synonymous added',syns,stop,missense,c;



for tx in TXS.keys(): print tx,TXS[tx]


#################################################################################################
## run tabix to get mutation data from EXAC VCF file
#tabix -p vcf /media/drive2/Variant-call-datasets/65Kexomes-DGM.v0.1/ExAC.r0.1.sites.vep.vcf.gz  20:47989510-48099180 > KCNB1_exac.tabix

tabix_command = ["/home/vbansal/Public/tools/htslib-master/tabix"];
tabix_command.append("/media/drive2/Variant-call-datasets/65Kexomes-DGM.v0.1/ExAC.r0.1.sites.vep.vcf.gz");
tabix_command.append(chrom + ':' + `firstpos` + '-' + `lastpos`); 
#tabix_command.append("> gene_EXAC.data");
#print tabix_command;
f = open('gene_EXAC.data','wb');
subprocess.call(tabix_command,stdout=f);
f.close();


# 20      47989752        .       G       A,C     17618.44        PASS    AC=3,1;AC_AFR=0,1;AC_AMR=0,0;AC_Adj=3,1;
EXAC_mutations = {}; 
File= open('gene_EXAC.data','r');
for line in File:
	var = line.strip().split('\t');
	chrom = var[0]; position = var[1]; ref = var[3]; alleles = var[4].split(',');
	info = var[7].split(';'); 
	AC = info[0].split('=')[1].split(','); 
	for i in xrange(len(alleles)):
		alt = alleles[i]
		try: 
			m = mutations[(var[0],int(var[1]),var[3],alt)]; 
			m[8]= 1; m[10] = int(AC[i])
			#print chrom,position,ref,alt,m;
		except KeyError: 
			pass; 
			#print 'notfound';
	#print var[0:5];


#################################################################################################
score = 0.87932;
score = 1.0;

for tx in TXS.iterkeys(): 
	sums = [0.0,0.0,0.0,0.0]; counts = [0,0,0,0];
	varlist1 = []; varlist2 = [];
	for mut in mutations.iterkeys():
		#print mutations[mut]
		if mutations[mut][8] == 1 and tx in mutations[mut][6]: 
			try: 
				if float(mutations[mut][3]) > score-epsilon: 
					sums[0] += mutations[mut][7]; counts[0] +=1; 
					varlist1.append(mutations[mut][10]);
					print mutations[mut],mut
				else: 
					varlist2.append(mutations[mut][10]);
				sums[1] += mutations[mut][7]; counts[1] +=1; 
			except ValueError: pass; 
		elif tx in mutations[mut][6]: 
			try: 
				if float(mutations[mut][3]) > score-epsilon: sums[2] += mutations[mut][7]; counts[2] +=1;
				sums[3] += mutations[mut][7]; counts[3] +=1;
			except ValueError: pass; 

	print tx,sums,counts;
	print varlist1;
	print varlist2;
