#last modified jan 23 2014

#cat ESP6500SI.allchroms.vcf | awk '{ split($8,F,";"); split(F[5],G,"="); split(G[2],A,","); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\tAF="A[2]";"F[3]; }' > ESP6500SI.allchroms.vcf.AA_AF 
#cat ESP6500SI.allchroms.vcf | awk '{ split($8,F,";"); split(F[5],G,"="); split(G[2],A,","); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\tAF="A[2]/100";"F[3]; }' > ESP6500SI.allchroms.vcf.AA_AF 
#python add_allelefreqs_VCF.py AF ../Variant-call-datasets/2000DanishExomes/2000exomes.AF.vcf phase1-190pools/combined.vcf.annotated

# grep VP= plist.515.ASN.CRISP.012514 | cut -f1-8 | grep -v LowDepth | awk '{ split($8,G,";"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"G[8]";"$8; }' > /projects/stsi/vbansal/T2D-pooledseq-july2012/AlleleFreq-data/plist.515.ASN.CRISP.012514.VCF

## script that takes a file with allele freqs and a VCF file and outputs all variants that overlap the regions in the bed file 

python add_allelefreqs_VCF.py AlleleFreq-data/2000exomes.Danish.AF.vcf phase1-190pools/combined.vcf.annotated > Danish.2000 &
python add_allelefreqs_VCF.py AlleleFreq-data/ESP6500SI.allchroms.vcf.AA_AF phase1-190pools/combined.vcf.annotated > AA.2250 &
python add_allelefreqs_VCF.py AlleleFreq-data/ESP6500SI.allchroms.vcf.EA_AF phase1-190pools/combined.vcf.annotated > EA.4400 &

python add_allelefreqs_VCF.py AlleleFreq-data/plist.505.EUR.CRISP.012814.VCF phase1-190pools/combined.vcf.annotated > EUR.1KG.505 &
python add_allelefreqs_VCF.py AlleleFreq-data/plist.507.AFR.CRISP.012814.VCF phase1-190pools/combined.vcf.annotated > AFR.1KG.507 &
python add_allelefreqs_VCF.py AlleleFreq-data/plist.515.ASN.CRISP.012514.VCF phase1-190pools/combined.vcf.annotated > ASN.1KG.515 &
paste EUR.1KG.505 ASN.1KG.515 AFR.1KG.507 EA.4400 AA.2250 Danish.2000 > Allpops.codingvars


## for Danish file, the two allele were switched so that the first allele was reference...

~                                                                                                             
