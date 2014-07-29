#!/usr/bin/python
import os, glob,sys, subprocess,re,random,math, gzip

##### take Ali's annotation phased file and generate haplotypes for each gene and calculate # of reads supporting each haplotype

def generate_gene_haplotypes(pileupfile):

        genes = {};  # genes is a hashtable indexed by gene-ids that stores the list of SNPs that overlap this gene 

        snptable = {};  # list of het-SNPs and indexes into the genelist....

        File = open(pileupfile,'r');
        for line in File:
                variant = line.split('\t');
                genelist = variant[16].split('///'); gtlist = variant[17].split('///');
                chrom = variant[1]; position = variant[3]; phased = variant[6]; a1 = variant[7]; a2 = variant[8];


                #if len(genelist) != len(gtlist): print genelist,gtlist,line;
                for i in xrange(len(genelist)):
                        if gtlist[i] == 'Noncoding_RNA': continue; # do we need this filter ??? 
                        try:
                                gene = genes[genelist[i]][0];
                                if gene[-1][1] != int(position): gene.append([chrom,int(position),a1,a2,phased,0,0]);

                        except KeyError:
                                genes[genelist[i]] = [[[chrom,int(position),a1,a2,phased,0,0]],0,0]; # reads supporting the two alleles

#               print gtlist,chrom,position,phased,a1,a2;

                genelist.sort();
                snptable[chrom + '_' + position] = [genelist[0]]; lastadded = genelist[0];
                for i in xrange(1,len(genelist)):
                        if genelist[i] != lastadded: snptable[chrom + '_' + position].append(genelist[i]); lastadded = genelist[i];

        File.close();
        return [genes,snptable];



###### add code to phase 'HET' SNPs prior to counting number of reads from rnaseq-hairs
###### if a read from RNA-seq data covers multiple SNPs, we change the 'HET' to RNA-PHASED with the two alleles appropriately phased....

##### use information about which strand a gene is on for read counting from RNA-seq hairs 

def calculate_haplotype_counts(rnaseqfile,genes,snptable):

        geneCU = {};
        File= open(rnaseqfile,'r');
        for line in File:
                read = line.strip().split();    chrom = 'chr' + read[4];

                ### set the gene update counter for each gene that is spanned by any variant in this read to 0 
                for i in xrange(5,len(read)): # loop over all the variants spanned by a read
                        try:
                                position = int(read[i].split('_')[0]); base = read[i].split('_')[1]; var = chrom + '_' + `position`;
                                genelist = snptable[var];
                                for gene in genelist: geneCU[gene] = 0;
                        except KeyError: print var,'not found'; pass;


                for i in xrange(5,len(read)): # loop over all the variants spanned by a read
                        try:
                                position = int(read[i].split('_')[0]); base = read[i].split('_')[1]; var = chrom + '_' + `position`;
                                genelist = snptable[var];
                                for i  in xrange(len(genelist)):  # loop over all genes covered by these variants 
                                        gene = genelist[i];
                                        try:
                                                genehap = genes[gene];
                                                for vars in genehap[0]:
                                                        if position == vars[1] and base == vars[2]:
                                                                vars[5] += 1;
                                                                if vars[4] != 'HET' and geneCU[gene] ==0: genehap[1] += 1; geneCU[gene] =1;
                                                                break;
                                                        elif position == vars[1] and base == vars[3]:
                                                                vars[6] += 1;
                                                                if vars[4] != 'HET' and geneCU[gene] == 0: genehap[2] += 1; geneCU[gene] = 1;
                                                                break;
                                                #print var,base,gene,genehap;
                                        except KeyError: pass;

                        except KeyError: print var,'not found'; pass;




        for gene,varlist in genes.iteritems():
                print "'"+gene+"'",varlist[1],varlist[2],
                for var in varlist[0]: print var[0]+':'+`var[1]` + ':'+ var[2] +':' + var[3] + ':' + var[4] + ':' + `var[5]` + ':' + `var[6]`,
                print;

        """
        for snp, genelist in snptable.iteritems():
                print snp,genelist;
        """


#####################################################################################################

# python create_gene_haplotypes.py CGI-data/phased_variant_annotations.ref.sort.Exon_3UTR_hets PC01_PolyA12.rnaseq.hairs.sorted.new > PC01_PolyA12.rnaseq.haplotypecounts 

[genes,snptable] = generate_gene_haplotypes(sys.argv[1]);
calculate_haplotype_counts(sys.argv[2],genes,snptable);
