#!/bin/bash

input_fastq_dir=<>
output_dir=<>
for mate1 in ${input_fastq_dir}/*_R1.fq
do
 mate2="${mate1/_R1.fq/_R2.fq}"
 samplename=${mate1%_R1.fq}

 echo Handling sample $samplename
 run-t1k -1 $mate1 -2 $mate2 --preset hla -f hlaidx/hlaidx_rna_seq.fa –od $output_dir -t 20

 echo "DONE!!"
done
