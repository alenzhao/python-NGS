#!/usr/bin/python
import os, glob,sys, subprocess,re,math,random 

## this code needs to be fixed to be used as a function call 

"""
1. read gene file from dbSNF and store list of all possible non-syn mutations in the gene... what about synonymous mutations ??
 python pvalue.py /media/drive2/Variant-call-datasets/dbNSFPv2.8.baylor/KCNB1.data ~/Public/tools/reference-genomes/NCBI37/human_g1k_v37.fasta
"""

from read_fasta_bychrom import make_fasta_index, read_chromosome


if len(sys.argv) < 3: print >>sys.stderr, "python program.py gene_name/file_name reference.fasta"; sys.exit()

## fasta_file 
offset_index = make_fasta_index(sys.argv[2]); current_chrom = '-';


if '/' in sys.argv[1]: dbSNFP_file = sys.argv[1]; 
else:
	gene_name = sys.argv[1]; 
	dbSNFP_file = gene_name + '.dbSNFP.data'; print >>sys.stderr, 'path',dbSNFP_file
	if not os.path.isfile(dbSNFP_file):
		# find chromosome for the gene 
		refseqfile = open('/home/vbansal/CODE/JOINTCODE-coral/PYTHON-scripts/VCF-programs/refseq.genes.ncbi37','r');
		for line in refseqfile: 
			tx = line.split(); 
			#print tx[12]
			if tx[12] == gene_name: chromosome = tx[2]; print >>sys.stderr, chromosome, gene_name; break; 

		## read the chromosome file for that gene 
		dirpath = "/media/drive2/Variant-call-datasets/dbNSFPv2.8.baylor/dbNSFP2.8_variant." + chromosome; 
		f = open(dbSNFP_file,'a');
		subprocess.call(["grep","#chr", dirpath],stdout=f);
		f.close();

		f = open(dbSNFP_file,'a');
		#print dirpath,gene_name 
		#subprocess.call(tabix_command,stdout=f);
		subprocess.call(["grep",gene_name,dirpath],stdout=f);
		f.close();


################################ read gene file from dbNSFP ######################################
mutations = {}; TXS = {}; firstpos= -1; lastpos = -1; mlist = {}; prevpos = 0; currentpos = -1;
syns = 0; missense = 0; stop = 0; c = 0; 
File = open(dbSNFP_file,'r'); 
for line in File: 
	var = line.split('\t');
	if line[0] == '#' or (var[4] != 'X' and var[5] != 'X' and var[4] != '.' and var[5] != '.'): 
		print var[0] + ':' + var[1] + ':' + var[2] + ':' + var[3],var[4] + '->' + var[5],var[22],var[23],var[29],var[30],var[31],
		print 'SIFT' + ':' + var[26] + ':' + var[28],
		print 'MTAST' + ':' + var[38] + ':' + var[39] + ':' + var[40],
		print 'GERP' + ':' + var[59] + ':' + var[60] + ':' + var[61]
	#print var[0],var[1],var[2],var[3],var[4],var[5],var[10],var[22],var[23],var[24],var[25],var[29],var[30],var[31]
	if line[0] == '#': continue;
	continue;

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

