#!/bin/bash

input_fastq_dir=~/fastq_files

for mate1 in ${input_fastq_dir}/*_R1.fq
do
 mate2="${mate1/_R1.fq/_R2.fq}"
 samplename=${mate1%_R1.fq}

 echo First alignment to whole genome
 hisat2 -p 20 -x ./genome --no-spliced-alignment -X 1000 -1 $mate1 -2 $mate2 -S ${samplename}.sam

 echo converting sam file to bam file

 samtools view -Sb ${samplename}.sam > ${samplename}_unsorted.bam
 echo removing sam file
 rm ${samplename}.sam


 echo sorting bam file by read_names
 samtools sort -n ${samplename}_unsorted.bam -o ${samplename}_sorted.bam -@ 20 #(sorts by read names)

 echo removing the previous temporary unsorted file

 rm ${samplename}_unsorted.bam

 echo filtering sorted bam file by proper read_pairs

 samtools view -b -f 2 ${samplename}_sorted.bam > ${samplename}_filtered.bam

 echo removing sorted bam file

 rm ${samplename}_sorted.bam

 echo converting filtered bam file into filtered fq files

 bedtools bamtofastq -i ${samplename}_filtered.bam -fq ${samplename}_filtered_R1.fq -fq2 ${samplename}_filtered_R2.fq

 echo genotyping HLA filtered fq files with HISATGENOTYPE

 python3 hisatgenotype --base hla --suffix fq -p 7 --pp 6 -1 ${samplename}_filtered_R1.fq -2 ${samplename}_filtered_R2.fq -v
 echo "DONE!"
done
