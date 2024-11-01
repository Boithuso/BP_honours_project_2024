#!/bin/bash
#GATK DepthOfCoverage module requires read groups or else its "WellformedReadFilter" filters
#all reads without RGs
input_bam_file_dir=<>
for bam_file in $input_bam_file_dir
do
 output_bam="${bamfile%.bam}_w_RG.bam"
 echo Adding read groups to $bam_file
 java -jar picard.jar AddOrReplaceReadGroups \
 I=$bamfile \
 O=$output_bam \
 RGID=RG1 \
 RGLB=Lib1 \
 RGPL=ILLUMINA \
 RGPU=PU1 \
 RGSM=Sample1
 
 echo Indexing $ouput_bam
 samtools index $output_bam
 echo "DONE!!"
done
