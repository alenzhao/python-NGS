#!/usr/bin/python
import os, glob,sys, subprocess,re,random,math, gzip

##### take a VCF file and refseq gene info and annotate the indels and SNPs
#### refseq file coordinates are 0-based while VCF is 1-based !!!
#### refseq file is stranded (genes on - strand go from right to left)

def read_transcripts(genefile):
	File = open(genefile,'r');
	TXlist = []; transcripts =0; chromosomeIndex = {};

	for line in File:
		if line[0] == '#': continue;
		transcript = line.strip().split('\t');
		chrom = transcript[2]; name = transcript[1]; strand = transcript[3]; 
		if 'hap' in chrom or 'random' in chrom: continue; 
		txstart = int(transcript[4]); txend = int(transcript[5]); cdstart = int(transcript[6]); cdend = int(transcript[7]);
		exons = int(transcript[8]);
		Exonstarts = transcript[9].rstrip(',').split(','); Exonends = transcript[10].rstrip(',').split(',');
		score = transcript[11]; hugoname = transcript[12]; 
		frames = transcript[15];
#		print hugoname,chrom,txstart,txend;
	
		TXlist.append([chrom,txstart,txend,cdstart,cdend,hugoname,name,exons,Exonstarts,Exonends,strand]); transcripts +=1;

	File.close();
	TXlist.sort();
	for i in xrange(transcripts): 
		if TXlist[i][0] not in chromosomeIndex: chromosomeIndex[TXlist[i][0]] = i; 
	print >>sys.stderr, "read",transcripts,'transcripts from file:',genefile;
	return [TXlist,transcripts,chromosomeIndex];

	for tx in TXlist: print tx;

def compare_variant_transcript(variant,transcript,annolist):
	## need to take strand information into account for saying Exon1 vs exon_last 
	if variant[1] >= transcript[1] and variant[1] < transcript[3]: 
		if transcript[10] == '+': annolist.append(['5UTR']);
		if transcript[10] == '-': annolist.append(['3UTR']);
	elif variant[1] > transcript[4] and variant[1] <= transcript[2]: 
		if transcript[10] == '+': annolist.append(['3UTR']);
		if transcript[10] == '-': annolist.append(['5UTR']);
	else: 
		flag =0;
		for i in xrange(transcript[7]):
			if variant[1] >= int(transcript[8][i]) and variant[1] <= int(transcript[9][i]): 
				if transcript[10] == '+': annolist.append(['Exon'+`(i+1)` + '/' + `transcript[7]`]);
				if transcript[10] == '-': annolist.append(['Exon'+`(transcript[7]-i)` + '/' + `transcript[7]`]);
				flag = 1; 
				break;
		if flag ==0:
			#for i in xrange(transcript[7]):
			annolist.append(['Intron']);
			
	#print transcript,

##INFO=<ID=GL,Number=.,Type=String,Description="geneList">	
####INFO=<ID=GM,Number=.,Type=String,Description="accession">
##FILTER=<ID=INDEL5,Description="Nearby 1000 Genomes Pilot Indels within 5bp">
##is it in last exon or first exon, print this information 
## DBSNP id 
## use SIFT/PolyPhen database 
## Exon/Intron/Upstream/Downstream etc 
## distance to end of beginning of transcript (for indels)
## distance from coding start site -2433  +433

def read_VCFfile(vcffile,TXlist,transcripts,chromosomeIndex):
	if vcffile != '-' and vcffile != 'sys.stdin': File = open(vcffile,'r');
	else: File = sys.stdin; 
	VCFlist = []; index = 0; prevchrom = "";
	
	for line in File: 
		if line[0] == '#': print line.strip(); continue;
		var = line.strip().split();
		#print len(var);
		chrom1 = var[0]; position = int(var[1]); ref = var[3]; alt = var[4]; score = var[5]; filter = var[6]; INFO = var[7]; 
		variant = [chrom1,position,ref,alt]; annolist = []; ovtxlist = []; txs =0;
		## increase index until we hit the right chromosome 
		if not chrom1.startswith('chr'): chrom = 'chr' + chrom1;
		else: chrom = chrom1;
		if prevchrom != chrom: 
			print >>sys.stderr, "processing variants on chrom",chrom;
			try: index= chromosomeIndex[chrom]; 
			except KeyError: 
				prevchrom = chrom;
				print line.strip(),'NOMATCH';
				continue; 

		s = index; 
		### variant position is less than the transcription start/end site of the first gene in sorted list 
		if position >= TXlist[index][2]: 
			while position >= TXlist[index][2] and index < transcripts and TXlist[index][0] == chrom: 
				index +=1; 
				if index >= transcripts: break;
			s +=1; 
		if index >= transcripts: print line.strip(); continue;
		
		if position < TXlist[index][1] and TXlist[index][10] == '+': 
			#ovtxlist.append(index); annolist.append(['Upstream']); txs +=1;
			pass;
		elif position < TXlist[index][1] and TXlist[index][10] == '-': 
			#ovtxlist.append(index); annolist.append(['Downstream']); txs +=1;
			pass;

		while position >= TXlist[index][1] and position <= TXlist[index][2] and index < transcripts and TXlist[index][0] == chrom: 
			compare_variant_transcript(variant,TXlist[index],annolist);
			ovtxlist.append(index); 
			index +=1; 
			if index >= transcripts: break;
			txs +=1;
		
		if txs ==0: 
			print line.strip();	
			#print >>sys.stderr, 'PROBLEM';
		else: 
			genelist = 'GL=' + ','.join([TXlist[ovtxlist[i]][5] for i in xrange(txs)]);
			NMlist = 'NM=' + ','.join([TXlist[ovtxlist[i]][6] for i in xrange(txs)]);
			ANlist = 'ANNO=' + ','.join([annolist[i][0] for i in xrange(txs)]);
			EXlist = 'EC=' + ','.join([`TXlist[ovtxlist[i]][7]` for i in xrange(txs)]);
			for i in xrange(7): print '%s\t' %(var[i]),
			print '%s\t' %(genelist + ';' + NMlist + ';' + ANlist + ';' + var[7]),
			for i in xrange(8,len(var)-1): print '%s\t' %(var[i]),
			print var[len(var)-1];
			#for i in xrange(txs): print TXlist[ovtxlist[i]][5]+':'+TXlist[ovtxlist[i]][6]+':'+`TXlist[ovtxlist[i]][7]`+':'+annolist[i][0],
		
		index = s;
		prevchrom = chrom;	
	
	if vcffile != '-' and vcffile != 'sys.stdin': File.close();
		
		#print index,TXlist[index];
	
		
		
		
	
[TXlist,transcripts,chromosomeIndex] = read_transcripts(sys.argv[1]);
VCFlist = read_VCFfile(sys.argv[2],TXlist,transcripts,chromosomeIndex); 


