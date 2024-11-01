!/bin/bash

#the input directory must contain bam files that are indexed and have there read groups
#within the gatk directory my have a list "input_bams.list" with the bam files

python3 gatk DepthOfCoverage -R genome.fa -O kir_DOC -I input_bams.list \
-DF NotSecondaryAlignmentReadFilter \
-gene-list kir_genes_and_gene_predictions.refseq \
-L KV575254.1 -L KV575246.1 -L KV575256.1 -L KV575253.1 -L KV575252.1 -L KV575255.1 -L KV575259.1 \
-L KI270917.1 -L KI270918.1 -L KI270919.1 -L KV575247.1 -L KV575248.1 -L KV575250.1 -L KV575249.1 \
-L KV575257.1 -L -L KI270920.1 -L KI270921.1 -L KI270922.1 -L KI270923.1 -L KI270929.1 -L KI270930.1 \
-L KI270931.1 -L KI270932.1 -L -L KI270933.1 -L KI270882.1 -L KI270883.1 -L KI270884.1 -L KI270885.1 \
-L KI270886.1 -L KI270887.1 -L KI270888.1 -L KV575258.1 -L KV575251.1 -L KV575260.1 -L KI270889.1 \
-L KI270890.1 -L GL000209.2 -L KI270891.1 -L KI270914.1 -L KI270915.1 -L KI270916.1 -L GL949746.1 \
-L GL949747.2 -L GL949748.2 -L GL949749.2 -L GL949750.2 -L GL949751.2 -L GL949752.1 -L GL949753.2 \
-L GL383573.1 -L -L GL383574.1 -L GL383575.2 -L KI270866.1 -L GL383576.1 -L KI270867.1 -L KI270865.1 \
-L KI270938.1 -L KI270868.1 -L KI270868.1 \ 
-L chr19:55200000-55390000
