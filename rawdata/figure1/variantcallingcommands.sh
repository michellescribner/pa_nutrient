
# Analysis of P. aeruginosa populations propagated in glucose, lactate, and amino acids for twelve days with and without biofilm selection

module purge
module load trimmomatic/trimmomatic-0.36
module load breseq/breseq-0.35.0 

# Ancestral clone
breseq -r /home/mrs186/pa14_nodrug/finalbreseq/breseqv35/Pseudomonas_aeruginosa_UCBPP-PA14_109.gbk /home/mrs186/pa14_nodrug/SRA/anc_forward_paired.fq.gz /home/mrs186/pa14_nodrug/SRA/anc_forward_unpaired.fq.gz /home/mrs186/pa14_nodrug/SRA/anc_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/SRA/anc_reverse_unpaired.fq.gz  -o /home/mrs186/pa14_nodrug/finalbreseq/breseqv35/anc -p --polymorphism-minimum-variant-coverage-each-strand 3 --consensus-minimum-variant-coverage-each-strand 3 -j 8

# Planktonic and biofilm populations 4 and 5
for i in B4 B5 P4 P5 ; do trimmomatic PE -phred33 -threads 4 /home/mrs186/pa14_nodrug/SRA/"$i"_*R1_001.fastq.gz /home/mrs186/pa14_nodrug/SRA/"$i"_*R2_001.fastq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_forward_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_forward_unpaired.fq.gz  /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/opt/trimmomatic/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70

breseq -r /home/mrs186/pa14_nodrug/finalbreseq/breseqv35/Pseudomonas_aeruginosa_UCBPP-PA14_109.gbk /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_forward_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_forward_unpaired.fq.gz  /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_reverse_unpaired.fq.gz -o /home/mrs186/pa14_nodrug/finalbreseq/breseqv35/"$i" -p --polymorphism-minimum-variant-coverage-each-strand 3 --consensus-minimum-variant-coverage-each-strand 3 -j 8 ; done

# Planktonic and biofilm populations 1-3
for i in 04 05 06 33 34 35 ; do trimmomatic PE -phred33 -threads 4 /home/mrs186/tobm9/reads_041818/041818_"$i"/*R1_001.fastq.gz /home/mrs186/tobm9/reads_041818/041818_"$i"/*R2_001.fastq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_forward_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_forward_unpaired.fq.gz  /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/opt/trimmomatic/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70

breseq -r /home/mrs186/pa14_nodrug/finalbreseq/breseqv35/Pseudomonas_aeruginosa_UCBPP-PA14_109.gbk /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_forward_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_forward_unpaired.fq.gz  /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_reverse_paired.fq.gz /home/mrs186/pa14_nodrug/finalbreseq/trimmed/"$i"_reverse_unpaired.fq.gz -o /home/mrs186/pa14_nodrug/finalbreseq/breseqv35/"$i" -p --polymorphism-minimum-variant-coverage-each-strand 3 --consensus-minimum-variant-coverage-each-strand 3 -j 8 ; done

# Parse outputs to xlsx file
/home/mrs186/scripts/BreseqCatEdited.py -p -d /home/mrs186/pa14_nodrug/finalbreseq/breseq35_3/

# copy outputs to computer
for i in B4 B5 P4 P5 04 05 06 33 34 35 anc ; do mkdir /Users/mrs/Documents/test/"$i"
scp -r mrs186@beagle.mmg.pitt.edu://home/mrs186/pa14_nodrug/finalbreseq/breseqv35/"$i"/output /Users/mrs/Documents/test/"$i" ; done

