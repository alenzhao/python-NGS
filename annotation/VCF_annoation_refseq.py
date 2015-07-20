#!/usr/bin/python
import os, glob,sys, subprocess,re,random,math, gzip

### author: Vikas Bansal
##### take a VCF file and refseq gene info and annotate the indels and SNPs ###


#### refseq file coordinates are 0-based while VCF is 1-based, refseq file is stranded (genes on - strand go from right to left)

### SEC61A exon from 127771678 to 127771745 (1 based) represented as 127771677-127771745 (last base is not exonic in 0-base) 
### python VCF_annoation_refseq.py refseq.genes.ncbi37 /media/drive2/T2D-pooledseq-july2012/newannotation-111914/phase2.vcf.8cols.annotated | cut -f1-9 > temp.out

UPSTREAM_WINDOW =2000; 

def calculate_AAlength(transcript):

	AAlength = 0.0; 
	for i in xrange(transcript[7]): 
		if int(transcript[9][i]) < transcript[3] or int(transcript[8][i]) > transcript[4] : AAlength += 0 
		elif int(transcript[8][i]) < transcript[3]: AAlength += int(transcript[9][i]) - transcript[3]; ## UTR
		elif int(transcript[9][i]) > transcript[4]: AAlength += transcript[4] - int(transcript[8][i]); ## UTR
		else: AAlength += int(transcript[9][i]) -int(transcript[8][i]);
	AAlength /= 3; 
	#print transcript[8],transcript[9],transcript[1:5],
	return AAlength;


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
		"""
                print transcript[1],hugoname,chrom,txstart,txend,cdstart,cdend,len(ExonStarts),
                for i in xrange(len(ExonStarts)):
                        s = int(ExonStarts[i]); e = int(ExonEnds[i]); 
                        #if s < cdstart: s = cdstart; 
                        #if e > cdend: e = cdend; 
                        print `s` + '-' + `e` + ':' + `(e-s+1)`,
                #print ExonStarts,ExonEnds,len(ExonStarts),len(ExonEnds);
                print;
                """
		TXlist.append([chrom,txstart,txend,cdstart,cdend,hugoname,name,exons,Exonstarts,Exonends,strand,0]); 
		AAlen = calculate_AAlength(TXlist[-1]); TXlist[-1][-1] = AAlen;
		transcripts +=1;

	File.close();
	TXlist.sort();
	for i in xrange(transcripts): 
		if TXlist[i][0] not in chromosomeIndex: chromosomeIndex[TXlist[i][0]] = i; 
	print >>sys.stderr, "read",transcripts,'transcripts from file:',genefile;
	return [TXlist,transcripts,chromosomeIndex];

	for tx in TXlist: print tx;

## compare variant to the gene 
def compare_variant_transcript(variant,transcript,annolist):
	## need to take strand information into account for saying Exon1 vs exon_last 
	position = variant[1]; lastposition = variant[4]; exons = transcript[7];
	FC = int(transcript[8][0]); ## first coding base of transcript
	if transcript[3] > FC: FC = transcript[3]; 
	LC = int(transcript[9][exons-1]);
	if transcript[4] < LC: LC = transcript[4]; 

	if (position >= transcript[1]) and lastposition < transcript[3]: ## between the start of transcript and first coding position....
		if transcript[10] == '+': annolist.append(['UTR5:-'+`(position-transcript[1])`]);
		if transcript[10] == '-': annolist.append(['UTR3:+' + `(transcript[3]-position)`]);
	elif position > transcript[4] and lastposition <= transcript[2]: 
		if transcript[10] == '+': annolist.append(['UTR3:+' + `(position-transcript[4])`]);
		if transcript[10] == '-': annolist.append(['UTR5:-' + `(transcript[2]-position)`]);
	else: 
		flag =0; 

		## transcript[1] is 5'UTR start and [2] is 3'UTR end | [3] is coding start and [4] is coding end...

		for i in xrange(exons): ## iterate over exons in the transcript...
			if position == int(transcript[9][i]) and transcript[10] == '+': print 'last base of exon+',position;
			if position == int(transcript[8][i])+1 and transcript[10] == '-': print 'last base of exon-',position;
			if (position <= int(transcript[8][i]) and lastposition > int(transcript[8][i]) and lastposition <= int(transcript[9][i]) ) or (position > int(transcript[8][i]) and position <= int(transcript[9][i]) and lastposition > int(transcript[9][i])): 
				interval = ',' + transcript[8][i] + ':' + transcript[9][i] + ',junction'; 
				print >>sys.stderr, 'variant',transcript[0],position,lastposition,variant[2],variant[3],'overlaps exon-intron junction/splicing for exon:',transcript[8][i],transcript[9][i];
				if transcript[10] == '+': annolist.append(['Exon'+`(i+1)` + '/' + `transcript[7]` + interval]);
				if transcript[10] == '-': annolist.append(['Exon'+`(transcript[7]-i)` + '/' + `transcript[7]` + interval]);
				flag = 1; 
				break;
			elif position > int(transcript[8][i]) and position <= int(transcript[9][i]): 
				
				interval = `(position-int(transcript[8][i])-1)` + '|.|' + `(int(transcript[9][i])-position)`; 
				if transcript[10] == '+': annolist.append(['Exon'+`(i+1)` + '/' + `transcript[7]` + ';POS=' + interval]);
				if transcript[10] == '-': annolist.append(['Exon'+`(transcript[7]-i)` + '/' + `transcript[7]` + ';POS=' + interval]);
				flag = 1; 
				break;
			elif position <= int(transcript[8][i]) and lastposition >= int(transcript[9][i]): 
				interval = ',' + transcript[8][i] + ':' + transcript[9][i] + ',span'; 
				print >>sys.stderr, 'variant',transcript[0],position,lastposition,variant[2],variant[3],'spans entire exon:',transcript[8][i],transcript[9][i];
				if transcript[10] == '+': annolist.append(['Exon'+`(i+1)` + '/' + `transcript[7]`+ interval]);
				if transcript[10] == '-': annolist.append(['Exon'+`(transcript[7]-i)` + '/' + `transcript[7]` + interval]);
				flag = 1; 
				break;
		if flag ==0:
			for i in xrange(exons-1): ## iterate over exons in the transcript...
				if position > int(transcript[9][i]) and position <= int(transcript[8][i+1]): 
					distance_5p = position-int(transcript[9][i]);
					distance_3p = int(transcript[8][i+1])-lastposition+1;
					interval = ';POS=' + `distance_5p` + '|.|' + `distance_3p`; 
					if transcript[10] == '+' and distance_5p < 8: annolist.append(['SpliceVar5:+' + `distance_5p`+ ':exon' + `(i+1)` + '/' + `transcript[7]`]);
					elif transcript[10] == '+' and distance_3p < 8: annolist.append(['SpliceVar3:-' + `distance_3p`+ ':exon' + `(i+2)` + '/' + `transcript[7]`]);
					if transcript[10] == '-' and distance_5p < 8: annolist.append(['SpliceVar3:-' + `distance_5p`+ ':exon' + `(exons-i+1)` + '/' + `transcript[7]`]);
					elif transcript[10] == '-' and distance_3p < 8: annolist.append(['SpliceVar5:+' + `distance_3p`+ ':exon' + `(exons-i)` + '/' + `transcript[7]`]);
					elif transcript[10] == '+' and (distance_3p <= 50 or distance_5p <= 50): annolist.append(['Intron'+`(i+1)` + '/' + `transcript[7]` + interval]);
					elif transcript[10] == '+' : annolist.append(['Intron'+`(i+1)` + '/' + `transcript[7]`]);
					elif transcript[10] == '-' and (distance_3p <= 50 or distance_5p <= 50): annolist.append(['Intron'+`(exons-i)` + '/' + `transcript[7]` + interval]);
					elif transcript[10] == '-': annolist.append(['Intron'+`(exons-i)` + '/' + `transcript[7]`]);

					print 'Intronic variant',position,distance_5p,distance_3p,transcript[10];
					flag = 1; 
					break;
			if flag ==0: 
				print >>sys.stderr, "variant not annotated as intronic",position,transcript[8][0],transcript[8][exons-1],transcript[10],transcript[1:5]
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
		chrom1 = var[0]; position = int(var[1]); lastposition = position;
		ref = var[3]; alt = var[4]; score = var[5]; filter = var[6]; INFO = var[7]; 
		lastposition += len(ref)-1;  # for big deletion 
		variant = [chrom1,position,ref,alt,lastposition]; annolist = []; ovtxlist = []; txs =0;
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
		### variant position is greater than the transcription start/end site of the first gene in sorted list 

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

		## indel is between gene start and end
		while (position >= TXlist[index][1] or lastposition >= TXlist[index][1]) and position <= TXlist[index][2] and index < transcripts and TXlist[index][0] == chrom: 
			compare_variant_transcript(variant,TXlist[index],annolist);
			ovtxlist.append(index); 
			index +=1; 
			if index >= transcripts: break;
			txs +=1;
		
		if txs ==0: 
			print line.strip();	
			#print >>sys.stderr, 'PROBLEM';
		else: 
			alleles = var[4].split(','); FS = [];
			for a in alleles: 
				len_var = len(var[3])-len(a); 
				if len_var < 0: len_var *= -1; 
				if len_var ==0: FS.append('SNV');
				elif len_var%3 ==0: FS.append('3n-indel');
				else: FS.append('frameshift')
			genelist = 'GL=' + ','.join([TXlist[ovtxlist[i]][5] for i in xrange(txs)]);
			NMlist = 'NM=' + ','.join([TXlist[ovtxlist[i]][6] for i in xrange(txs)]);
			ANlist = 'ANNO=' + ','.join([annolist[i][0] for i in xrange(txs)]);
			EXlist = 'EC=' + ','.join([`TXlist[ovtxlist[i]][7]` for i in xrange(txs)]);
			Lengths = 'AL=' + ','.join([`TXlist[ovtxlist[i]][11]` for i in xrange(txs)]);
			FSlist = 'FS=' + ','.join([FS[i] for i in xrange(len(alleles))]); 
			for i in xrange(7): print '%s\t' %(var[i]),
			if 'Exon' in ANlist: print '%s\t' %(FSlist + ';' + genelist + ';' + NMlist + ';' + ANlist + ';' + Lengths + ';\t' + var[7]), ## frameshift additional annotation
			else: print '%s\t' %(genelist + ';' + NMlist + ';' + ANlist  + ';\t' + var[7]),
			print ANlist + ';' + Lengths + ';\t' + var[7],
			for i in xrange(8,len(var)-1): print '%s\t' %(var[i]),
			print var[len(var)-1];
			#for i in xrange(txs): print TXlist[ovtxlist[i]][5]+':'+TXlist[ovtxlist[i]][6]+':'+`TXlist[ovtxlist[i]][7]`+':'+annolist[i][0],
		
		index = s;
		prevchrom = chrom;	
	
	if vcffile != '-' and vcffile != 'sys.stdin': File.close();
		
		#print index,TXlist[index];
	
		
		
		
	
[TXlist,transcripts,chromosomeIndex] = read_transcripts(sys.argv[1]);
VCFlist = read_VCFfile(sys.argv[2],TXlist,transcripts,chromosomeIndex); 


