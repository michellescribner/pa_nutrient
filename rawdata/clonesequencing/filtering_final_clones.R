### Filtering for P. aeruginosa clones chosen for sequencing from populations propagated in CF media for 12 days

library(tidyverse)

setwd("/Users/mrs/Documents/pa14_nodrug/submit/rawdata/clonesequencing")

#save SNP tab as csv
#Convert to noutf characters 
system("iconv -c -f utf-8 -t ascii//TRANSLIT Breseq_Output_snps.csv > breseq_output_snps_noutf.csv")
snps <- read.csv("breseq_output_snps_noutf.csv",header=TRUE)
snps$Gene <- gsub(" <- ", "", snps$Gene)
snps$Gene <- gsub(" -> ", "", snps$Gene)
snps$Gene <- gsub(" <-", "", snps$Gene)
snps$Gene <- gsub(" ->", "", snps$Gene)
nrow(snps)#4732

#include sample names
samplekey <- read.csv("samplekey_pa14nodrug_clones.csv",header=TRUE)
snps <- (merge(snps,samplekey,by="Sample"))

ancestor_snps <- subset(snps,Sample == "anc")
## Remove all ancestral mutations based on the Position column
snps_noref <- snps[ !(snps$Position %in% ancestor_snps$Position), ]
nrow(snps_noref)#156

snps_p <- pivot_wider(snps_noref, id_cols = c(Annotation, Gene, Description), names_from = Sample, values_from = Evidence, values_fill = NA)

### Manually filter mutations
# Certain mutations occurred in most populations, but weren't called in the ancestral clone's analysis. 
# However, based on their trajectories they likley were in the ancestral clone but merely didn't meet the criteria to be called 
snps_noref <- subset(snps_noref, snps_noref$Gene != "PA14_13130/PA14_13140") 
snps_noref <- subset(snps_noref, snps_noref$Gene != "PA14_16820/PA14_16830")
snps_noref <- subset(snps_noref, snps_noref$Gene != "gcd/PA14_34990")
snps_noref <- subset(snps_noref, snps_noref$Gene != "pvcA/ansA")
snps_noref <- subset(snps_noref, snps_noref$Gene != "fabI/ppiD")
snps_noref <- subset(snps_noref, snps_noref$Gene != "rpsF/PA14_65190") 
nrow(snps_noref)

write.csv(snps_noref,file=("snps_noref.csv"))

######################

### New Junction Evidence
system("iconv -c -f utf-8 -t ascii//TRANSLIT Breseq_Output_nje.csv > breseq_output_nje_noutf.csv")
nje <- read.csv("breseq_output_nje_noutf.csv",header=TRUE)

#nje splits mutations into two rows (one row for each side of the junction), must combine into one row for further analysis
evens <- seq(from = 2, to = nrow(nje), by = 2)
nje_even <- nje[evens, ]
cols_even <- paste(colnames(nje_even), "2", sep = ".")
colnames(nje_even) <- cols_even
nje_odd <- nje[-evens, ]
nje <- cbind(nje_odd, nje_even)

#include sample names
nje <- (merge(nje,samplekey,by.y="Sample", by.x="Sample"))

nje$Position <- gsub('"', "", as.character(nje$Position))
nje$Position <- gsub("\\^A", "", as.character(nje$Position))
nje$Position.2 <- gsub('"', "", as.character(nje$Position.2))

# remove mutations in the ancestor
nje <- subset(nje, nje$Gene != "PA14_35700/PA14_35710") #this is likely circularization of an IS element
# not supported by 3 reads on each strand
nje <- subset(nje, nje$Gene != "nrdA") 
nje <- subset(nje, nje$Gene != "PA14_56910") 

write.csv(nje, file="nje.csv")

