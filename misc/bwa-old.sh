#!/bin/bash
#PBS -A stsi-group
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l mem=24000M
cd /projects/stsi/PROJECTS/hla-bwa/;
genome=/home/vbansal-scripps/GENOMES/NCBI37/human_g1k_v37.fasta
#genome=/home/vbansal-scripps/GENOMES/1000genomes-huref/human_b36_male.fa
bwa=/home/vbansal-scripps/bin/bwa-0.5.9/bwa
DATA=/projects/stsi/PROJECTS/HLA/raw

#$bwa index -a bwtsw $genome;

for sample in NA07348_L6 NA07348_L7 NA07348_L8 #NA10835_L6 NA10835_L7 NA10835_L8 NA10835_L6 NA10835_L7 NA10835_L8 NA12750_L6 NA12750_L7 NA12750_L8
do
echo $sample
ls -l $DATA/$sample\_1.fastq
$bwa aln -t 8 -L  -I  $genome $DATA/$sample\_1.fastq > $sample\_1.sai
$bwa aln -t 8 -L  -I  $genome $DATA/$sample\_3.fastq > $sample\_3.sai
#$bwa aln -n 8 -o 4 -e 2 -t 8 -L -k 3 -I  $genome $DATA/$sample\_1.fastq > $sample\_1.sai
#$bwa aln -n 8 -o 4 -e 2 -t 8 -L -k 3 -I  $genome $DATA/$sample\_3.fastq > $sample\_3.sai
$bwa sampe $genome $sample\_1.sai $sample\_3.sai $DATA/$sample\_1.fastq $DATA/$sample\_3.fastq > $sample.sam
/home/vbansal-scripps/bin/samtools-0.1.7a/samtools import $genome.fai $sample.sam $sample.bam
#/home/vbansal-scripps/bin/samtools-0.1.7a/samtools sort $sample.bam $sample.sorted
rm -f $sample\_1.sai $sample\_3.sai $sample.sam
done

