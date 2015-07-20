import sys,os
from optparse import OptionParser

from read_muttable import get_mutation_rates

## translate gene into protein sequence... gene as list of exons, genome as text file, translation table 3bases-> amino acid 
## take a transcript or gene-model (list of intervals), chrom name, output sequences corresponding to the gene (appropriate strand)...
## given frameshift indel, find new protein sequence and location of Stop codon 
## determine if by skipping an exon or even two exons, we get a functional transcript in frame... 
## for all possible locations of frameshift indels (1-10 bp), determine impact... 
## use knowledge of alternate splicing to determine all possibilities... 

## frameshift + coding SNP -> avoid stop codon... interesting haplotype effects... 

from read_fasta_bychrom import make_fasta_index, read_chromosome
from codontable import gencode, AAnames, revcomp

## IMPORTANT: exon starts, ends are 0-based and [a,b) interval format 
## for last exon, calculate extra translation in different frames... 
VNTRtable = {}; 

## also allow some mismatches in VNTR tract.. likelihood model = copy of VNTR + mismatch 
def calculate_rightshift(chromosome,start,end,maxindel,chrom_id,gene_id):
	VNTRlist = []; prevlength = 1; prevshift = -1; prevstart = 0; prevend= 0;
	## for each base in exon, examine all deletions of length 1-K for VNTR/MS
	for position in xrange(start,end+1):
		for lengthvar in xrange(1,maxindel+1):

			leftpos = position; rightpos = position+lengthvar; rightshift = 0;
			while chromosome[0][leftpos] == chromosome[0][rightpos]:  leftpos +=1; rightpos +=1; rightshift +=1; 
			if rightshift >= lengthvar*2 and rightpos > prevend: 
				try: 
					info = VNTRtable[(chrom_id,position,lengthvar,rightpos)];
				except KeyError: 
					VNTRtable[(chrom_id,position,lengthvar,rightpos)] = 1; 
					if lengthvar > 6: 
						print 'VNTR',gene_id,chrom_id,position,'unit-length',lengthvar,chromosome[0][position:position+lengthvar],rightshift+lengthvar,
						#print chromosome[0][position:rightpos],
						print 'COPIES',rightshift/lengthvar+1,chrom_id+':'+`position`+ '-' + `(rightpos+1)`; 
					elif rightshift/lengthvar >= 4: 
						print 'repeat_TR',gene_id,chrom_id,position,'unit-length',lengthvar,chromosome[0][position:position+lengthvar],rightshift+lengthvar,
						print chromosome[0][position-2:rightpos+2],
						print 'COPIES',rightshift/lengthvar+1,chrom_id+':'+`position`+ '-' + `(position+1)`; 
				
				prevlength = lengthvar; prevshift = rightshift; prevstart = position; prevend = rightpos; 




## for a given transcript, output the coding sequence and translated protein sequence 
def output_exon_sequences(chromosome,CDS_start,CDS_end,exon_starts,exon_ends,exons,strand,TRANSLATE_FLAG,chrom,gene_id,PFLAG):
	## chromosome is simply the two tuple (sequence,length)
	## strand is for outputting + or -
	leftoverseq = ''; codingseqlist = []; Tlist = [];
	for i in xrange(exons):
		first = exon_starts[i]; last = exon_ends[i]; 
		if first < CDS_start: first= CDS_start; 
		if last > CDS_end: last = CDS_end; 

		if last < CDS_start or first > CDS_end: ## non-coding exon
			exonseq = ''; ls = 0; leftoverseq = '';
			if PFLAG: print 'Exon',i+1,strand,exon_starts[i],exon_ends[i],ls,'leftover:',ls%3, chromosome[0][first:last],'NON_CODING_EXON';
		else: ## coding exon 

			if strand == '+': exonseq = leftoverseq + chromosome[0][first:last]; ls = len(exonseq); exonseq1 = chromosome[0][first:last]; 
			else: 
				beans = chromosome[0][first:last]; ls = len(beans);
				exonseq1 = ''.join([revcomp.get(beans[ls-j-1],'X') for j in xrange(ls)]);
				exonseq = leftoverseq + ''.join([revcomp.get(beans[ls-j-1],'X') for j in xrange(ls)]);
				ls += len(leftoverseq); 

			translated = ''.join([ gencode.get(exonseq[j:j+3],'-') for j in xrange(0,ls-2,3)]); SC = 0; 
			if PFLAG: print 'Exon',i+1,strand,exon_starts[i],exon_ends[i],ls,'L%3='+ `ls%3`,leftoverseq, exonseq1#,translated;


			if TRANSLATE_FLAG =="1":
				translated1 = ''.join([ gencode.get(exonseq[j:j+3],'-') for j in xrange(1,ls-2,3)]); SC1 = 0;
				for j in xrange(len(translated1)): 
					if translated1[j] == '*': SC1 +=1; 
				translated2 = ''.join([ gencode.get(exonseq[j:j+3],'-') for j in xrange(2,ls-2,3)]); SC2 = 0;
				for j in xrange(len(translated2)): 
					if translated2[j] == '*': SC2 +=1; 

				if PFLAG: 
					print '#PROT#',translated;
					print '#PROT#',translated1,
					if SC1 ==0: print '#NSC#', 
					print '+1';
					print '#PROT#',translated2,
					if SC2 ==0: print '#NSC#',
					print '+2';
				#for j in xrange(0,ls-2,3): print exonseq[j:j+3],gencode.get(exonseq[j:j+3],'-'),j;
				calculate_rightshift(chromosome,first,last,50,chrom,gene_id);
			
			if ls%3 ==0: leftoverseq = '';
			elif ls%3 ==1: leftoverseq = exonseq[-1]; 
			elif ls%3 ==2: leftoverseq = exonseq[ls-2:ls];
			
			## add flanking bases for exon, useful for mutation rate estimation 
			if strand == '+': codingseqlist.append([chromosome[0][first:last],first,last,last-first,chromosome[0][first-1],chromosome[0][last],chrom,strand]); 
			else: 
				beans = chromosome[0][first:last]; ls = len(beans);
				rbeans = ''.join([revcomp.get(beans[ls-j-1],'X') for j in xrange(ls)]);
				codingseqlist.append([rbeans,last-1,first-1,last-first,revcomp.get(chromosome[0][last],'X'),revcomp.get(chromosome[0][first-1],'X'),chrom,strand]);
			Tlist.append(translated); 

	len_protein = -1;
	proteinseq = ''.join(Tlist); 
	for i in xrange(len(proteinseq)): 
		if proteinseq[i] == '*' and len_protein == -1: len_protein = i; 

	#if PFLAG: print >>sys.stderr, 'FULL_PROTEIN:',len_protein,'\t',''.join(Tlist);
	#for p in Tlist: len_protein += len(p); 

	return [codingseqlist,proteinseq,len_protein]; 

## output set of all mutations of a given type (stop missense, synonymous) 
## codingseqlist is an array of 4 tuples (codingseq (exon), start,end, strand), start and end correspond to + strand, codingseq is correct strand
def output_all_mutations(codingseqlist,mutation_type,outfilename):

	### read file with information about probability of mutation A[B]C -> A[B*]C where B->B*
	mut_table = get_mutation_rates();
	if outfilename != "None": outFile = open(outfilename,'w');
	else: outFile = sys.stdout;

	print '\n ###########function to analyze coding seq for stop gain mutations#########\n ';

	codingseq = ''.join([codingseqlist[i][0] for i in xrange(len(codingseqlist))]);
	proteinseq = ''.join([gencode.get(codingseq[j:j+3],'-') for j in xrange(0,len(codingseq)-2,3)])
	l = len(codingseq);
	#print 'FULL_PROTEIN',proteinseq,'len_coding_seq',l,l%3;

	# need to output exon number for each mutation based on exon lengths 
	#for c in codingseqlist: print c[1],c[2],c[3];
	curr_exon = 0; t = 0;
	for i in xrange(0,l):
		if t ==0: trimer = codingseqlist[curr_exon][4] + codingseq[i] + codingseq[i+1]; 
		elif t == codingseqlist[curr_exon][3]-1: trimer = codingseq[i-1] + codingseq[i] + codingseqlist[curr_exon][5]; 
		else: trimer = codingseq[i-1:i+2]; 

		frame = i/3;
		codon = codingseq[frame*3:frame*3+3]; AA = gencode.get(codon,'-')
		for m in ['A','C','G','T']: 
			if m == codingseq[i]: continue; 
			if i == frame*3: newcodon = m + codingseq[frame*3+1:frame*3+3]; 
			elif i == frame*3+1: newcodon = codingseq[frame*3] + m + codingseq[frame*3+2]; 
			else: newcodon = codingseq[frame*3] + codingseq[frame*3+1] + m;
			newAA = gencode.get(newcodon,'-');
			pflag = 0;
			if codingseqlist[curr_exon][2] > codingseqlist[curr_exon][1]: position =  codingseqlist[curr_exon][1]+t
			else: position = codingseqlist[curr_exon][1]-t
			if (mutation_type == "stop" or mutation_type == "nonsense") and newAA == "*": pflag = 1; 	
			elif (mutation_type == "missense" or 'nonsyn' in mutation_type) and newAA != AA and newAA != '*': pflag =1;
			elif (mutation_type == "silent" or mutation_type == 'synonymous') and newAA == AA and newAA != '*': pflag =1;
			elif (mutation_type == "all"): pflag = 1; 
			if newAA == AA: TYPE = 'silent'; 
			elif newAA == "*" and AA != '*': TYPE = 'stopgain'; 
			elif AA == "*" and newAA != AA: TYPE = 'stoplost';
			elif newAA != AA: TYPE = 'missense'; 

			if codingseqlist[curr_exon][7] == '+': ref = codingseq[i]; alt = m; 
			else: ref = revcomp[codingseq[i]]; alt = revcomp[m]; 

			if pflag ==1: print >>outFile, 'Exon',curr_exon+1,i+1,frame+1,codon,newcodon,codingseq[i], m,AA,newAA,'trimer',trimer,m,t,codingseqlist[curr_exon][3],mut_table[(trimer,m)],codingseqlist[curr_exon][6] + ':' + `position+1` + ':' + ref + ':' + alt,'p.' + AAnames[AA]+`(frame+1)`+AAnames[newAA],TYPE;
		t +=1; 
		if t >= codingseqlist[curr_exon][3]: t -= codingseqlist[curr_exon][3]; curr_exon +=1; 

	if outfilename != "None": outFile.close();

	

### new function added 05/06/2015 
def process_transcript(Tx,line): # Tx is line split into columns where line read from refseq.ncbi37 file 

	strand = Tx[3]; Txstart = int(Tx[4]); Txend = int(Tx[5]);  CDS_start = int(Tx[6]); CDS_end= int(Tx[7]); exons = int(Tx[8]);
	exon_starts = [int(i) for i in Tx[9].rstrip(',').split(',')]; 
	exon_ends = [int(i) for i in Tx[10].rstrip(',').split(',')]; 
	if strand == '-': exon_starts.reverse(); exon_ends.reverse();

	print '######',Tx[1],Tx[12],chrom,strand,exons,Txstart,Txend,CDS_start,CDS_end,'#####';
	print line,
	[codingseqlist,proteinseq,len_protein] = output_exon_sequences(CHROM,CDS_start,CDS_end,exon_starts,exon_ends,exons,strand,options.translate,chrom,Tx[12],1); 
	print;
	if options.muttype != "None": output_all_mutations(codingseqlist,options.muttype,options.outfilename);


def find_canonical_transcript(transcripts):

	nt = len(transcripts);
#	for i in xrange(1,nt):
#		if transcripts[i][0] > transcripts[longest_tx][0]: longest_tx = i;

	canonical = -1; 
	for t in xrange(nt):
		Tx = transcripts[t][1]; 
		strand = Tx[3]; Txstart = int(Tx[4]); Txend = int(Tx[5]);  CDS_start = int(Tx[6]); CDS_end= int(Tx[7]); exons = int(Tx[8]);
		exon_starts = [int(i) for i in Tx[9].rstrip(',').split(',')]; 
		exon_ends = [int(i) for i in Tx[10].rstrip(',').split(',')]; 
		if strand == '-': exon_starts.reverse(); exon_ends.reverse();
		[codingseqlist,proteinseq,len_protein] = output_exon_sequences(CHROM,CDS_start,CDS_end,exon_starts,exon_ends,exons,strand,0,chrom,Tx[12],0); 
		if len_protein > 0: print transcripts[t][1][1],len_protein,transcripts[t][0],Txstart,Txend,CDS_start,CDS_end
		transcripts[t].append(len_protein);
		if canonical < 0 or len_protein > transcripts[canonical][3] or (len_protein == transcripts[canonical][3] and exons > int(transcripts[canonical][1][8])) or  (len_protein == transcripts[canonical][3] and exons == int(transcripts[canonical][1][8]) and abs(Txstart-Txend) > abs(int(transcripts[canonical][1][5])-int(transcripts[canonical][1][4]))): canonical = t; 
	print >>sys.stderr, 'canonical transcript',transcripts[canonical][1][1]
	return canonical	
	

#######################################################################################################################

## two arguments, arg1 is reference.fasta arg2 is transcript model file refseq.genes.table

parser = OptionParser();
parser.add_option("--ref",dest="reference",help="reference fasta file",default="None");
parser.add_option("--model",dest="genemodelFile",help="refseq transcript file",default="None");
parser.add_option("--gene",dest="gene",help="gene to parse",default="None");
parser.add_option("--translate",dest="translate",help="if set to 1, translate the coding sequence per exon",default="0");
parser.add_option("--mutation_type","--type",dest="muttype",help="",default="None");
parser.add_option("--out",dest="outfilename",help="",default="None");
(options,args) = parser.parse_args();


gene_name = options.gene; 
if options.reference == 'None' or options.genemodelFile == 'None': 
	print >>sys.stderr," need two arguments, --ref=reference.fasta file AND  --model=refseq.transcript file ";
	print >>sys.stderr," optional third argument is gene_name (MLL2, GCK, etc) ";
	sys.exit();
else: 
	print >>sys.stderr,'reference:',options.reference;
	print >>sys.stderr,'translate genes',options.translate; 
	#sys.exit()

offset_index = make_fasta_index(options.reference); 
current_chrom = '-';

################# read refseq transcript file and process each transcript #####################

File = open(options.genemodelFile); 
if not File: print >>sys.stderr, 'genemodel file does not exist'; sys.exit();
transcripts = [];
for line in File: 
	if line[0] == '#': continue;
	Tx = line.strip().split();
	if gene_name != "None" and Tx[12] != gene_name and Tx[1] != gene_name and Tx[2] != gene_name: continue; 

	chrom = Tx[2].strip('chr');  

	if chrom != current_chrom: ## read in new chromosome 
		if chrom in offset_index: CHROM  = read_chromosome(options.reference,offset_index,chrom); 
		else: continue;
		current_chrom = chrom;
	
	if gene_name == "None": process_transcript(Tx,line); 
	else: transcripts.append([int(Tx[8]),Tx,line]);


## choose longest transcript if there are multiple transcripts, if we want another transcript, we can specify the RefSeqID/ENSEMBL-ID
if gene_name != "None" and len(transcripts) > 0:
	nt = len(transcripts); longest_tx = 0;
	if nt > 1: 
		longest_tx = find_canonical_transcript(transcripts)
		print >>sys.stderr, "\nThe Gene",gene_name,"has", nt,"transcripts:",
		for i in xrange(nt): print >>sys.stderr, transcripts[i][1][1] + ':' + `transcripts[i][0]`,
		print >>sys.stderr, "\nmutations will be output for only the longest transcript",transcripts[longest_tx][1][1]," !! Please specify transcript ID instead of gene name for analyzing a particular transcript\n\n";

	process_transcript(transcripts[longest_tx][1],transcripts[longest_tx][2]);





######################################################################


