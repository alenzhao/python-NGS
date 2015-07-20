##script for running annovar and generating excel file with pvalues for association and annotation for all variants ##
## last modified dec 16 2012 ## 
## author Vikas Bansal ## 

## had to modify the summarize_annovar code to change "," to "|": new code in file summarize_annovar_newdelim.pl

## make VCF file with 8 cols for running annovar 
cut -f1-8 136genes.variants.oct25.vcf.FILTERED > 136genes.variants.oct25.vcf.FILTERED.annovar

## make variant format file for annovar 
perl convert2annovar.pl --format vcf4 ../T2D-pooledseq-july2012/VariantCalls-oct082012/OCT25-vcf/136genes.variants.oct25.vcf.FILTERED.annovar --includeinfo > ../T2D-pooledseq-july2012/VariantCalls-oct082012/OCT25-vcf/ANN

## run ANnotation pipeline using annovar script 
perl summarize_annovar_changedelim.pl --ver1000g 1000g2012feb --buildver hg19 --verdbsnp 135 ../T2D-pooledseq-july2012/VariantCalls-oct082012/OCT25-vcf/annotation/ANN humandb-hg19/

## manually add header to vcf FILTERED FILE before next step for association analysis
python ../../addannotations-new.py annotation/ANN.genome_summary.csv 136genes.variants.oct25.vcf.FILTERED >  136genes.variants.oct25.vcf.FILTERED.annotated

/home/vbansal-scripps/python/bin/python2.7 /home/vbansal-scripps/PROGRAMS/JOINTCODE-coral/PYTHON-scripts/NGS-library-python/genetests-sept72012.py 136genes.variants.oct25.vcf.FILTERED.annotated ../../Phenotypes-pools/192pools.phenotypes 40 > pvalues

cat pvalues| awk 'BEGIN { ORS = "\t"; }{ split($5,A,"/"); f = A[1]/A[2]; print $13"\t"$14"\t"$16"\t"$17"\t"$1"\t"$6"\t"$8"\t"$9"\t"$10"\t"f"\t"$2"\t"$3"\t"$4"\t"$5"\t"$15; printf "\n"; }' > pvalues.csv

cat pvalues | awk 'BEGIN { ORS = "\t"; }{ split($5,A,"/"); f = A[1]/A[2]; print $13"\t"$14"\t"$16"\t"$17"\t"$1"\t"$6"\t"$8"\t"$9"\t"$10"\t"f"\t"$2"\t"$3"\t"$4"\t"$5"\t"$12"\t"$15; printf "\n"; }' > pvalues.biandtriallelic.csv

cat pvalues | awk 'BEGIN { ORS = "\t"; }{ split($15,AN,";"); split($5,A,"/"); f = A[1]/A[2]; print $13"\t"$14"\t"$16"\t"$17"\t"$1"\t"$6"\t"$8"\t"$9"\t"$10"\t"f"\t"$2"\t"$3"\t"$4"\t"$5"\t"$12"\t"AN[1]"\t"AN[2]"\t"AN[3]"\t"AN[4]"\t"AN[5]"\t"AN[6]; printf "\n"; }' > pvalues.biandtriallelic.6colsplit.csv
## additional column for tri-allelic variants 

## 133 variants with minor allele frequency >= 0.05 and not in DBSNP -> triallelic or bogus variants (indels in homopolymer runs...) 

cat 136genes.variants.oct25.vcf.FILTERED.annotated | grep -v '^#' | awk 'BEGIN { ORS = "\t"; } {split($8,I,";");  for (i=1;i<=6;i++) print $i; print I[2]";"I[5]; for (i=10;i<NF;i++) print $i; printf "%s\n",$NF; }' > 136genes.variants.oct25.vcf.FILTERED.annotated.csv

grep exonic 136genes.variants.oct25.vcf.FILTERED.annotated.csv > 136genes.variants.oct25.vcf.FILTERED.annotated.coding.csv

## remove X chromosome p-values from association 

