#!/bin/bash

#SBATCH --job-name=clones
#SBATCH -o clonebreseq.out
#SBATCH --mail-user=mrs186@pitt.edu
#SBATCH --mail-type=FAIL,END

# Trimming and variant calling for P. aeruginosa clones isolated from day 12 populations

module purge
module load trimmomatic/trimmomatic-0.36
module load breseq/breseq-0.35.0 

for i in 120 121 122 125 126 127 128 129 130 131 ; do trimmomatic PE -phred33 /home/mrs186/pa14_nodrug/seq_200128/reads/012820_$i/*R1_001.fastq.gz /home/mrs186/pa14_nodrug/seq_200128/reads/012820_$i/*R2_001.fastq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_forward_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_forward_unpaired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/opt/trimmomatic/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70; done

for i in 25 26 27 28 29 30 31 32 33 34 35 36 ; do trimmomatic PE -phred33 /home/mrs186/pa14_nodrug/seq_191012/reads/101219_$i/*R1_001.fastq.gz /home/mrs186/pa14_nodrug/seq_191012/reads/101219_$i/*R2_001.fastq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_forward_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_forward_unpaired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/opt/trimmomatic/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70; done

trimmomatic PE -phred33 /home/mrs186/pa14_nodrug/SRA/1021_1_S86_R1_001.fastq.gz /home/mrs186/pa14_nodrug/SRA/1021_1_S86_R2_001.fastq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/1021_1_forward_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/1021_1_forward_unpaired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/1021_1_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/1021_1_reverse_unpaired.fq.gz ILLUMINACLIP:/opt/trimmomatic/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70

for i in 1021_1 120 121 122 125 126 127 128 129 130 131 25 26 27 28 29 30 31 32 33 34 35 36 ; do breseq -r /home/mrs186/pa14_nodrug/finalbreseq/breseqv35/Pseudomonas_aeruginosa_UCBPP-PA14_109.gbk /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_forward_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_forward_unpaired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq_clones/trimmed/"$i"_reverse_unpaired.fq.gz -o /home/mrs186/pa14_nodrug/finalbreseq_clones/breseq/"$i" --consensus-minimum-variant-coverage-each-strand 3 -j 4; done

breseq -r /home/mrs186/pa14_nodrug/finalbreseq/breseqv35/Pseudomonas_aeruginosa_UCBPP-PA14_109.gbk /home/mrs186/pa14_nodrug/SRA/anc_forward_paired.fq.gz /home/mrs186/pa14_nodrug/SRA/anc_forward_unpaired.fq.gz /home/mrs186/pa14_nodrug/SRA/anc_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/SRA/anc_reverse_unpaired.fq.gz -o /home/mrs186/pa14_nodrug/finalbreseq_clones/breseq/anc --consensus-minimum-variant-coverage-each-strand 3 -j 4

module load miniconda
/home/mrs186/scripts/BreseqCat_clones.py -d /home/mrs186/pa14_nodrug/finalbreseq_clones/breseq
/home/mrs186/scripts/BreseqCatCoverage.py -d /home/mrs186/pa14_nodrug/finalbreseq_clones/breseq





