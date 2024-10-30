#!/bin/bash

#Setting up a virtual machine (VM) for analysis with HISAT-genotype, T1K, and Geny
#VM specs:
#region: us-central1(iowa)
# n2-highmem-32(32 vCPU, 16 core, 256 GB memory)
# - 32 threads 
# - boot disk 150 GB is sufficient if you plan on continually loading and off loading onto and off of the vm
#Access scopes: set access for each API: allow use of google storage buckets
#after activating the vm:


#installing essentials:

echo Installing necessary tools

sudo apt-get update
sudo apt-get install git
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install sysstat
sudo apt-get install samtools
sudo apt-get install bedtools
sudo apt install python3.11-venv
sudo apt-get install zlib1g-dev
sudo apt-get install python3-dev

echo --------------------------------------------------------------------------------------
#Hisatgenotype installation:
#HISAT-genotype unfortunatetly only be installed and run from the home directory, also running HISAT-genotype commands nessitates that you specify the absolute path to the home directory and not just use the tilde ~

echo "Building Hisat-genotype and exporting to PATH"

git clone --recurse-submodules https://github.com/DaehwanKimLab/hisat-genotype ~/hisatgenotype

cd ~/hisatgenotype/hisat2

make

cd ..

export PATH=~/hisatgenotype:~/hisatgenotype/hisat2:$PATH
export PYTHONPATH=~/hisatgenotype/hisatgenotype_modules:$PYTHONPATH


echo running hisatgenotype setup.sh to download pre-built references
bash setup.sh -r

echo downloading and indexing ensemble release 112 

cd indicies
rm genome.fa
wget ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
samtools faidx genome.fa

echo Building hla.graph references
python3 ../hisat2/hisat2-build -p 5 --snp hla.index.snp --haplotype hla.haplotype hla_backbone.fa hla.graph
cd ..

echo "Performing genotyping on HISAT-genotype test dataset"
## Performing this test not only not makes sure the tool is working but also allows the HLA database to be built ###

echo Downloading test dataset
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/hla/ILMN.tar.gz
tar xvzf ILMN.tar.gz

echo Genotyping ILMN/NA12892 test data with hisatgenotype
python3 hisatgenotype --base hla --locus-list A -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz -v

###Output should look like:

                        #1496 reads and 769 pairs are aligned
                                #1 A*02:01:01:01 (count: 419)
                                #2 A*02:251 (count: 407)
                                #3 A*02:610:02 (count: 405)
                                #4 A*02:524:02 (count: 403)
                                #5 A*02:650 (count: 403)
                                #6 A*02:01:123 (count: 402)
                                #7 A*02:647 (count: 401)
                                #8 A*02:01:126 (count: 400)
                                #9 A*02:562 (count: 400)
                                #10 A*02:645 (count: 400)


                                #1 ranked A*02:01:01:01 (abundance: 51.95%)
                                #2 ranked A*11:01:01:01 (abundance: 48.05%)

echo --------------------------------------------------------------------------------------

#For geny:
cd ~

git clone https://github.com/0xTCG/geny
cd geny
python3 -m venv geny_venv
source geny_venv/bin/activate
pip install numpy
pip install pysam
pip install psutil
pip install mappy
pip install pyyaml
pip install tqdm

deactivate


echo --------------------------------------------------------------------------------------
#For T1K:
cd ~

echo Installing T1K 

git clone https://github.com/mourisl/T1K.git
cd T1K

echo Building allele reference sequences from IPD-IMGT/HLA
perl t1k-build.pl -o hlaidx --download IPD-IMGT/HLA -d https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat -g gencode.gtf

echo Building allele reference sequences from IPD-KIR
perl t1k-build.pl -o kiridx --download IPD-KIR --partial-intron-noseq -d https://ftp.ebi.ac.uk/pub/databases/ipd/kir/kir.dat -g gencode.gtf

#echo Genotyping example dataset with T1K
#./run-t1k -f kiridx/kiridx_rna_seq.fa -1 example/example_1.fq -2 example/example_2.fq -t 8 -o T1K_example

echo "DONE!!"
cd
