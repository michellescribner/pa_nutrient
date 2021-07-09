### Filtering for Evolution of PA14 in M9 Media Using the Bead Model: no drug populations ###

library(tidyverse)
library(ggrepel)

setwd("/Users/mrs/Documents/pa14_nodrug/submit/rawdata/figure1")

# populations 1,2,3,4,and 5 were reanalyzed together using trimmomaticv36 and breseqv35 requiring 3 reads on each strand to support variants
# used PA14 genome 109 from winsor 2016
# the directories of Breseq data were parsed into an .xlsx file with tabs for SNPs,MC, and NJE, then copied to computer using
# /home/mrs186/scripts/BreseqCatEdited.py -d /home/mrs186/pa14_nodrug/finalbreseq/breseqv35

### SNPs
system("iconv -c -f utf-8 -t ascii//TRANSLIT Breseq_Output_snps.csv > Breseq_Output_snps_noutf.csv")
snps <- read.csv("Breseq_Output_snps_noutf.csv",header=TRUE)
snps$gene <- gsub(" <- ", "", snps$gene)
snps$gene <- gsub(" -> ", "", snps$gene)
snps$gene <- gsub(" <-", "", snps$gene)
snps$gene <- gsub(" ->", "", snps$gene)
snps$position <- gsub(":1", "", snps$position)
nrow(snps) #2835

colnames(snps) <- c("Sample", "Evidence", "Position", "Mutation", "Frequency", "Annotation", "Gene", "Description")
# include sample names
samplekey <- read.csv("samplekey_pa14_nodrug.csv",header=TRUE)
snps <- (merge(snps,samplekey,by="Sample"))
snps$Sample <- snps$Population

# remove % sign from frequency column
snps$Frequency <- gsub( "%", "", as.character(snps$Frequency))
# convert frequency values to numeric
snps$Frequency <- as.numeric(as.character(snps$Frequency))

### Remove SNPs already present in the ancestor
ancestor_snps <- subset(snps,Sample == "PA14 Ancestor")
snps_noref <- snps[ !(snps$Position %in% ancestor_snps$Position), ]
nrow(snps_noref)#138

### Remove manually filtered mutations
# multiple mutations within reads, poor mapping likely
snps_noref <- subset(snps_noref, snps_noref$Gene != "rpsF/PA14_65190") 
snps_noref <- subset(snps_noref, snps_noref$Gene != "PA14_61200")
# variant lacked support in ancestor but was present
snps_noref <- subset(snps_noref, snps_noref$Gene != "PA14_13130/PA14_13140") 
snps_noref <- subset(snps_noref, snps_noref$Gene != "gcd/PA14_34990")
snps_noref <- subset(snps_noref, snps_noref$Gene != "fabI/ppiD")
# unlikely trajectories, potentially in ancestor
snps_noref <- subset(snps_noref, snps_noref$Gene != "PA14_16820/PA14_16830")
snps_noref <- subset(snps_noref, snps_noref$Gene != "pvcA/ansA")
# junction not supported by 3 reads on each strand
snps_noref <- subset(snps_noref, snps_noref$Gene != "PA14_19940")
snps_noref <- subset(snps_noref, snps_noref$Gene != "PA14_06790")

snps_final <- separate(snps_noref, col = Population, sep = ", ", into = c("day", "pop", "b_p"), remove=FALSE)
# write.csv(snps_final,file=("snps_final.csv"))

# cast data frame - organizing with each mutation as the rows and the frequency of that mutation on a given day as the columns
snps_cast <- snps_final %>%
  pivot_wider(id_cols=c(Gene, Annotation, Mutation, Description, Position), names_from = Population, names_sort=TRUE, values_from = Frequency, values_fn = sum, values_fill = 0)
write.csv(snps_cast, file="snps_cast.csv")

snps_cast_gene <- snps_final %>%
  pivot_wider(id_cols=c(Gene, Description), names_from = Population, names_sort=TRUE, values_from = Frequency, values_fn = sum, values_fill = 0)
# write.csv(snps_cast_gene, file="snps_cast_gene.csv")

###############################################

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
nje <- (merge(nje,samplekey,by.y="Sample", by.x="sample"))

nje$freq <- gsub("%", "", as.character(nje$freq))
#convert frequency values to numeric
nje$freq <- as.numeric(as.character(nje$freq))

nje$position <- gsub('"', "", as.character(nje$position))
nje$position <- gsub("\\^A", "", as.character(nje$position))
nje$position.2 <- gsub('"', "", as.character(nje$position.2))

# write.csv(nje, file="nje.csv")

nje$info <-  paste(nje$position, nje$gene,nje$annotation, nje$product, nje$position.2, nje$gene.2,nje$annotation.2, nje$product.2, sep=":")
nje$info <- gsub('"', '', nje$info)

nje_cast <- nje %>%
  pivot_wider(id_cols=info, names_from = Population, names_sort=TRUE, values_from = freq, values_fn = sum, values_fill = 0)
nje_cast <- separate(nje_cast, col = info, sep = ":", into = c("Position1", "Gene1", "Annotation1", "Description1", "Position2", "Gene2", "Annotation2", "Description2"), extra="drop")
write.csv(nje_cast, file="nje_cast.csv")

########################################

### Figure 1 

# clean NJE to combine with SNPs
njemap <- nje %>% 
  select(Population, num, position, freq, annotation, gene, product)
colnames(njemap) <- c("Sample", "num", "Position", "Frequency", "Annotation", "Gene", "Description")
njemap$Position <- gsub("= ", "", as.character(njemap$Position))
njemap$Position <- gsub(" =", "", as.character(njemap$Position))

#remove mutations that were detected in the ancestor
ancestor_nje <- subset(njemap,Sample == "PA14 Ancestor")
njemap <- njemap[ !(njemap$Position %in% ancestor_nje$Position), ]

#remove mutations that represent prophage circularization because frequency is not useful and these may not reflect actual "mutations"
njemap <- subset(njemap, njemap$Position != "4345126")
njemap <- subset(njemap, njemap$Position != "4338302")

njemap$MutationType <- "New Junction"
snps_final$MutationType <- "SNP"
snps_final$MutationType[grep("\\*", snps_final$Annotation)] <- "Nonsense"
snps_final$MutationType[grep("coding", snps_final$Annotation)] <- "Small Indel"
snps_final$MutationType[grep("+AGCGGCGTCCCGCGAT", snps_final$Mutation)] <- "Small Indel"

snps_final <- snps_final %>%
  select(Sample, num, Position, Frequency, Annotation, Gene, Description, MutationType)
snps_final <- rbind(snps_final, njemap)

snps_final %>%
  select(Position, Annotation, Gene, Description, MutationType) %>% 
  distinct() %>%
  group_by(MutationType) %>%
  summarize(n = n())

snps_final$Function <- "other"
snps_final$Function[grep("lasR", snps_final$Gene)] <- "quorum sensing"
snps_final$Function[grep("4081494", snps_final$Position)] <- "quorum sensing"
snps_final$Function[grep("4071618", snps_final$Position)] <- "quorum sensing"
snps_final$Function[grep("morA", snps_final$Gene)] <- "cyclic-di-GMP"
snps_final$Function[grep("wsp", snps_final$Gene)] <- "cyclic-di-GMP"
snps_final$Function[grep("PA14_50060", snps_final$Gene)] <- "cyclic-di-GMP"
snps_final$Function[grep("pf5r", snps_final$Gene)] <- "Pf5 prophage"
snps_final$Function[grep("PA14_48810", snps_final$Gene)] <- "Pf5 prophage"
snps_final$Function[grep("pil", snps_final$Gene)] <- "type IV pili"
snps_final$Function[grep("PA14_64050", snps_final$Gene)] <- "cyclic-di-GMP"

snps_final$Position <- as.numeric(gsub(",", "", snps_final$Position))
snps_final$Sample <- gsub(", Day 12", "", snps_final$Sample)
snps_final$Sample <- gsub(", Population", "", snps_final$Sample)

c <- snps_final %>%
  mutate(Sample = fct_reorder(Sample, num)) %>%
  mutate(MutationType = fct_relevel(MutationType, "SNP", "Small Indel", "Nonsense")) %>%
  mutate(Function = fct_relevel(Function, "quorum sensing", "cyclic-di-GMP", "Pf5 prophage", "type IV pili", "other")) %>%
  ggplot(aes(y=Sample, x=Position)) + 
  geom_point(aes(color=Function, size=Frequency, shape=MutationType, stroke=2)) +
  scale_shape_manual(name = "Mutation Type", values=c(1, 2, 8, 0))+ 
  scale_color_manual(values = c("quorum sensing" = "#4393C3", "cyclic-di-GMP" = "#339966", "Pf5 prophage" = "#9933FF", "type IV pili" = "#404788FF", "other" = "dark gray")) +
  scale_size(name = "Allele Frequency", range = c(0.25,6), limits = c(0,100)) +
  scale_x_continuous(limits = c(0,6537648), breaks = c(2e6, 4e6, 6e6),label = c("2,000,000", "4,000,000", "6,000,000"), expand = c(0,0))+
  theme_minimal()+
  theme(legend.position="bottom", legend.box = "vertical", legend.title = element_text(size=12)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  theme(panel.grid.minor = element_blank(), legend.spacing.y = unit(-.1, "in")) +
  xlab("Genome Position (bp)") + ylab("Evolved Population")+
  guides(colour = guide_legend(override.aes = list(size=3)), shape = guide_legend(override.aes = list(size=3)))
#+ theme(legend.position = "none")
ggsave(filename= "Figure1.pdf", c, dpi=300, dev='pdf', width = 9, height = 5)

parallel <- snps_final %>% 
  group_by(Gene, Description) %>% 
  summarize(n()) 

#########################################

### Figure 2 ###

snps_final$exp <- "Laboratory"
snps_final$mut <- snps_final$Annotation
snps_final$mut <- gsub("\\s*\\([^\\)]+\\)","",as.character(snps_final$mut))

# Figure 2A

genes <- read.csv("/Users/mrs/Documents/pa14_nodrug/submit/rawdata/Figure2/Pseudomonas_aeruginosa_UCBPP-PA14_109_geneannotations.csv",header=TRUE)
genes <- genes %>% 
  filter(Start > 4070000 & Start < 4140000) 
genes %>% 
  ggplot(aes(x=Start, y=7, label = Gene.Name)) + 
  geom_rect(aes(xmin=Start, xmax=End, ymin=6, ymax=7), fill="gray")+
  scale_x_continuous(limits = c(4070000, 4140000), position = "top")+
  scale_y_continuous(limits = c(0, 7))+
  theme_minimal() +  xlab("Genome Position (bp)") +
  theme(legend.position = "none", axis.line = element_blank(), axis.text.y = element_blank(), 
        panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  geom_rect(aes(xmin=4081494, xmax=4130628, ymin=0.5, ymax=1), fill="#4393C3")+
  geom_rect(aes(xmin=4071618, xmax=4112326, ymin=1.5, ymax=2), fill="#4393C3")+
  geom_text_repel(aes(x=Start, y=6), direction = "x", angle = 90, max.overlaps = Inf, min.segment.length = 0, fontface = "italic", nudge_y = -1.5, size = 3.6) 
ggsave(filename= "Figure2a.1.pdf", dpi=300, dev="pdf", width = 9, height = 1.5, units="in")

snps_final %>% 
  filter(Position > 4050000 & Position < 4150000) %>% 
  ggplot(aes(y="This Study", x=Position, label = mut)) + 
  geom_vline(xintercept = 4085339)+
  geom_vline(xintercept = 4086058)+
  geom_point(aes(y = 1, x= Position, color=Function, shape=MutationType, stroke=2)) +
  scale_shape_manual(name = "Mutation Type", values=c("SNP" = 1, "Small Indel" = 2, "Nonsense" = 8, "New Junction" = 0)) + 
  scale_color_manual(values = c("quorum sensing" = "#4393C3")) +
  scale_x_continuous(limits = c(4085339, 4086058))+
  scale_y_continuous(limits = c(0, 2), breaks = 1)+
  annotate(geom = "rect", xmin = 4086058-(18*3), xmax = 4086058-(159*3), ymin = 0, ymax = Inf, fill = "navy blue", alpha = 0.1) +
  annotate(geom = "rect", xmin = 4086058-(176*3), xmax = 4086058-(231*3), ymin = 0, ymax = Inf, fill = "#4393C3", alpha = 0.15) +
  theme_minimal() +  xlab("Genome Position") +
  theme(legend.position = "none", panel.grid.minor = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=12))+
  geom_text_repel(aes(x=Position, y=1), direction = "x", angle = 60, nudge_y =0.05, hjust= "left", min.segment.length = Inf, force = 0.1, size =3.6) 
ggsave(filename= "Figure2a.2.pdf",dpi = 300, device="pdf", width = 8, height = 2)

# Figure 2B
morA <- read.csv("/Users/mrs/Documents/pa14_nodrug/submit/rawdata/Figure2/morA_snps_clinical.csv",header=TRUE)
morA$Function <- "cyclic-di-GMP"
morA$study_number <- as.character(morA$study_number)

morA$Domain <- "none"
morA$Domain[between(morA$aa_position, 305, 349)] <- "PAS"
morA$Domain[between(morA$aa_position, 608, 699)] <- "PAS"
morA$Domain[between(morA$aa_position, 729, 844)] <- "PAS"
morA$Domain[between(morA$aa_position, 865, 966)] <- "PAS"
morA$Domain[between(morA$aa_position, 978, 1139)] <- "DGC"
morA$Domain[between(morA$aa_position, 1159, 1394)] <- "PDE"

morA <- morA %>% 
  filter(synonymous == 0) %>%
  distinct(study_number, aa_position, aminoacid, Domain)
morA%>% 
  ggplot(aes(y=study_number, x=aa_position, label = aminoacid)) + 
  geom_point(aes(color=Domain, stroke=2), shape=1) +
  scale_color_manual(values = c("none" = "dark gray", "PAS" = "DCE319FF", "DGC" = "#73D055FF", "PDE" = "#3CBB75FF")) +
  scale_x_continuous(limits = c(0, 1416))+
  scale_y_discrete(limits = c("3", '2', "1", "NA"), labels = c("Clinical Isolates", "Laboratory Isolates", "This Study", "NA")) +
  theme_minimal() + theme(panel.grid.minor = element_blank(), text = element_text(size = 12))+ 
  xlab(element_blank()) + ylab(element_blank()) +
  geom_text_repel(direction = "x", angle = 60, nudge_y =0.05, hjust= "left", min.segment.length = Inf, force = 0.1, size =3.4) +
  theme(legend.position="bottom", legend.box = "vertical") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.position = "none")+
  annotate(geom = "rect", xmin = 305, xmax = 349, ymin = 0, ymax = Inf, fill = "#DCE319FF", alpha = 0.15) +
  annotate(geom = "rect", xmin = 608, xmax = 699, ymin = 0, ymax = Inf, fill = "#DCE319FF", alpha = 0.15) +
  annotate(geom = "rect", xmin = 729, xmax = 844, ymin = 0, ymax = Inf, fill = "#DCE319FF", alpha = 0.15) +
  annotate(geom = "rect", xmin = 865, xmax = 966, ymin = 0, ymax = Inf, fill = "#DCE319FF", alpha = 0.15) +
  annotate(geom = "rect", xmin = 978, xmax = 1139, ymin = 0, ymax = Inf, fill = "#73D055FF", alpha = 0.2) +
  annotate(geom = "rect", xmin = 1159, xmax = 1394, ymin = 0, ymax = Inf, fill = "#3CBB75FF", alpha = 0.2) 
ggsave(filename= "Figure2b.pdf", dpi=300, device = "pdf", width = 8, height = 4)

# Figures were then cropped and combined to produce final versions
# The snps_cast.csv and nje_cast.csv files were combined into Supplementary Data 1
