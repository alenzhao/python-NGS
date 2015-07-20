sample=$1;
cores=12
#/home/vbansal/bin/misc/sratoolkit.2.3.2-5-centos_linux64/bin/fastq-dump -A $sample --split-files /oasis/tscc/scratch/vbansal/NA12878-comparison/$sample.sra
/home/vbansal/bin/bwa-0.7.5a/bwa mem -M -t $cores -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" /projects/stsi/vbansal/T2D-pooledseq-july2012/target/ncbi37.fa $sample\_1.fastq.gz $sample\_2.fastq.gz  1> $sample.bwamem.sam 2> $sample.bwamem.log

ST=/home/vbansal/bin/samtools-0.1.19/samtools

$ST view -@ 4 -bS $sample.bwamem.sam > $sample.bwamem.bam;
$ST sort -@ 4  $sample.bwamem.bam $sample.bwamem.sorted;
java -Xmx2g -jar /projects/stsi/tools/picard/latest/MarkDuplicates.jar I=$sample.bwamem.sorted.bam O=$sample.bwamem.MD.bam M=$sample.picard.metrics AS=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/projects/stsi/vbansal/novocraft/temp CREATE_INDEX=true;
mv $sample.bwamem.MD.bai $sample.bwamem.MD.bam.bai






java -Xmx2g -jar /projects/stsi/tools/picard/latest/MergeSamFiles.jar I=N76629543_1.bwamem.sorted.bam I=N76629543_2.bwamem.sorted.bam O=N76629543.bwamem.merged.bam MSD=true &
java -Xmx2g -jar /projects/stsi/tools/picard/latest/MergeSamFiles.jar I=N990300_1.bwamem.sorted.bam I=N990300_2.bwamem.sorted.bam O=N990300.bwamem.merged.bam MSD=true &
java -Xmx2g -jar /projects/stsi/tools/picard/latest/MergeSamFiles.jar I=T76629543_1.bwamem.sorted.bam I=T76629543_2.bwamem.sorted.bam O=T76629543.bwamem.merged.bam MSD=true  &
java -Xmx2g -jar /projects/stsi/tools/picard/latest/MergeSamFiles.jar I=T990300_1.bwamem.sorted.bam I=T990300_2.bwamem.sorted.bam O=T990300.bwamem.merged.bam MSD=true  &



/home/vbansal/bin/misc/sratoolkit.2.1.7-centos_linux64/fastq-dump --split-files $sample.sra --gzip
