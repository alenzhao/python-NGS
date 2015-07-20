import sys

## postdoc test, write code to translate gene into protein sequence... gene as list of exons, genome as text file, translation table 3bases-> amino acid 

# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec12
#https://kaspermunch.wordpress.com/2013/11/19/finding-open-reading-frames/ 
#https://www.biostars.org/p/55851/
import urllib 
code = "Q7Z7W5"
#data = urllib.urlopen("http://www.uniprot.org/uniprot/" + code + ".txt").read()
#print data

#http://stackoverflow.com/questions/620367/python-how-to-jump-to-a-particular-line-in-a-huge-text-file
## read fasta file, build simple index, next time reuse it by writing to file, allows reading of one chromosome at a time... good enough for processing data... 
File= open(sys.argv[1],'r');
line_offset = []
offset = 0;
for line in File:
	if line[0] == '>': line_offset.append([offset,line.strip().split()[0]]); print line, 
	offset += len(line)
print >>sys.stderr, 'finished reading fasta file'
File.seek(16909); print 'line',File.readline(); 
File.seek(line_offset[20][0]); print 'line',File.readline(); 
while 1: 
	line = File.readline()
	if line[0] == '>': break; 
#print line_offset
#File.seek(16571,0);
File.close();
"""
"""

## take a transcript or gene-model (list of intervals), chrom name, output sequences corresponding to the gene (appropriate strand)...
## given frameshift indel, find new protein sequence and location of Stop codon 
## determine if by skipping an exon or even two exons, we get a functional transcript in frame... 
## for all possible locations of frameshift indels (1-10 bp), determine impact... 
## use knowledge of alternate splicing to determine all possibilities... 

beans = "TGACTGTGTTTCTGAACAATAAATGACTTAAACCAGGTATGGCTGCCGATGGTTATCTT"

revcomp = {'A':'T', 'T':'A', 'C':'G', 'G':'C','a':'T', 't':'A', 'c':'G', 'g':'C' } 

gencode = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
		'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def translate_frameshifted( sequence ):
	translate = ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])
	return translate
	    
print translate_frameshifted(beans[0:])    # first frame
print translate_frameshifted(beans[1:])    # second frame
print translate_frameshifted(beans[2:])    # third frame
print translate_frameshifted(beans[::-1][0:])    # first frame from end
print translate_frameshifted(beans[::-1][1:])    # second frame from end
print translate_frameshifted(beans[::-1][2:])    # third frame from end

l = len(beans);
rcbeans = ''.join([revcomp.get(beans[l-i-1],'X') for i in xrange(l)]);
print beans,'\n',beans[::-1],'\n',rcbeans;
