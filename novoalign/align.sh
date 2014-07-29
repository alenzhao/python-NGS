sample=$1;
FQ1=/projects/stsi3/PROJECTS/T2D-pooledseq/June5_2014/;
FQ2=/projects/stsi3/PROJECTS/T2D-pooledseq/May30_2014/;
NV=/projects/stsi/vbansal/novocraft/;
ST=/home/vbansal/bin/samtools-0.1.19/samtools
cores=8;
index=/projects/stsi/vbansal/novocraft/ncbi37.novoindex
outdir=/oasis/tscc/scratch/vbansal/T2D_phase3/$1;
mkdir -p $outdir;

#zcat $FQ1/Sample_$1/$1_*_L00[1-8]_R1_001.fastq.gz  > $outdir/$1_R1.fastq; 
#zcat $FQ1/Sample_$1/$1_*_L00[1-8]_R2_001.fastq.gz > $outdir/$1_R2.fastq;
#$NV/novoalign -c $cores --rOQ -r Random -d $index -H -f $outdir/$1_R1.fastq $outdir/$1_R2.fastq -o SAM -i PE 200,100 --softclip 40 -k -K $outdir/$1.calibration > $outdir/$1.novoalign.sam.1;

zcat $FQ2/Sample_$1/$1_*_L00[1-8]_R1_001.fastq.gz  > $outdir/$1_R1.fastq;
zcat $FQ2/Sample_$1/$1_*_L00[1-8]_R2_001.fastq.gz > $outdir/$1_R2.fastq;
$NV/novoalign -c $cores --rOQ -r Random -d $index -H -f $outdir/$1_R1.fastq $outdir/$1_R2.fastq -o SAM -i PE 200,100 --softclip 40 -k -K $outdir/$1.calibration > $outdir/$1.novoalign.sam.2;


#######################################################################################################

## align whole-genome sample

$NV/novoalign -c $cores --rOQ -r Random -d $index -H -f $sample\_1.fastq.gz $sample\_2.fastq.gz -o SAM "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA\tPI:300" --softclip 40  1> $sample.novoalign.sam 2> $sample.novoalign.log





########################################################################################################

ST=/home/vbansal/bin/samtools-0.1.19/samtools;
sample=$1;
outdir=/oasis/tscc/scratch/vbansal/T2D_phase3/$sample;

$ST view -bS $outdir/$sample.novoalign.sam.1 > $outdir/$sample.novoalign.bam.1;
$ST sort -m 3600M $outdir/$sample.novoalign.bam.1 $outdir/$sample.novoalign.sorted.1
$ST view -bS $outdir/$sample.novoalign.sam.2 > $outdir/$sample.novoalign.bam.2;
$ST sort -m 3600M $outdir/$sample.novoalign.bam.2 $outdir/$sample.novoalign.sorted.2
java -Xmx2g -jar /projects/stsi/tools/picard/latest/MergeSamFiles.jar I=$outdir/$sample.novoalign.sorted.1.bam I=$outdir/$sample.novoalign.sorted.2.bam O=$outdir/$sample.novoalign.merged.bam MSD=true
java -Xmx2g -jar /projects/stsi/tools/picard/latest/MarkDuplicates.jar I=$outdir/$sample.novoalign.merged.bam O=$outdir/$sample.novoalign.merged.MD.bam M=$outdir/$sample.picard.metrics AS=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/projects/stsi/vbansal/novocraft/temp CREATE_INDEX=true;
mv $outdir/$sample.novoalign.merged.MD.bai $outdir/$sample.novoalign.merged.MD.bam.bai


