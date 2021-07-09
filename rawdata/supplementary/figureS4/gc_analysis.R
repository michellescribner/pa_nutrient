### Analysis of growth curves of lasR and morA mutants in media supplemented with various carbon sources

# Three independent experiments were performed on 2.5.20, 2.12.20, and 7.24.20

library(readxl)
library(tidyverse)
library(growthcurver)
library(forcats)
library(multcompView)
library(car)

setwd("/Users/mrs/Documents/pa14_nodrug/submit/rawdata/supplementary/figureS4")

# read in data
data1 <- read_excel("gc_analysis.xlsx", sheet = 1)
colnames(data1) <- paste("all", colnames(data1), sep = "ppp")
data2 <- read_excel("gc_analysis.xlsx", sheet = 2)
colnames(data2) <- paste("glucose", colnames(data2), sep = "ppp")
data3 <- read_excel("gc_analysis.xlsx", sheet = 3)
colnames(data3) <- paste("amino acids", colnames(data3), sep = "ppp")
data4 <- read_excel("gc_analysis.xlsx", sheet = 4)
colnames(data4) <- paste("lactate", colnames(data4), sep = "ppp")
df <- as.data.frame(cbind(data1,data2,data3,data4))
df <- df %>%
  rename(time = allppptime)

# remove these rows because they contain outlier values for several samples
df <- df[-c(1,4),]
# this sample was an outlier thus was removed (high sigma value, >.1
df$'allpppmorA N1124Kppp7' <- NULL

### Growth curve analysis

num_analyses <- length(names(df)) - 1
d_gc <- data.frame(sample = character(num_analyses),
                     k = numeric(num_analyses),
                     n0  = numeric(num_analyses),
                     r = numeric(num_analyses),
                     t_mid = numeric(num_analyses),
                     t_gen = numeric(num_analyses),
                     auc_l = numeric(num_analyses),
                     auc_e = numeric(num_analyses),
                     sigma = numeric(num_analyses),
                     stringsAsFactors = FALSE)
d_fit <- vector("list", num_analyses)
names(d_fit) <- names(df[-1])
trim_at_time <- 24 # Change as necessary, would like to see the curves first.

n <- 1    # keeps track of the current row in the output data frame
for (col_name in names(df)) {
  if (col_name != "time") {
    # Create a temporary data frame that contains just the time and current col
    d_loop <- df[, c("time", col_name)]

    # Now, call Growthcurver to calculate the metrics using SummarizeGrowth
    gc_fit <- SummarizeGrowth(data_t = d_loop[, "time"], 
                              data_n = d_loop[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "min")
      
    # Now, add the metrics from this column to the next row (n) in the 
    # output data frame, and increment the row counter (n)
    d_gc$sample[n] <- col_name
    d_gc[n, 2:9] <- c(gc_fit$vals$k,
                      gc_fit$vals$n0,
                      gc_fit$vals$r,
                      gc_fit$vals$t_mid,
                      gc_fit$vals$t_gen,
                      gc_fit$vals$auc_l,
                      gc_fit$vals$auc_e,
                      gc_fit$vals$sigma)
    d_fit[col_name] <- gc_fit["data"]
    n <- n + 1
    }
  }

  
#look for outliers with bad fits (large sigma values)
gc_out <- as_tibble(d_gc)
# Plot a histogram of the sigma values in order to check for outliers
hist(gc_out$sigma, main = "Histogram of sigma values", xlab = "sigma")

# Show the top 5 samples with the largest sigma value 
# (with the worst model fit to the growth curve data)
gc_out %>% top_n(10, sigma) %>% arrange(desc(sigma))

out <- separate(gc_out, sample, into = c("media", "strain", "replicate"), sep = "ppp")

out$strain <- gsub("lasR.R216Q", "lasR R216Q", out$strain)
out$strain <- gsub("PA14_45920..PA14_46440.Δ49.133.b", "∆PA14_45920..PA14_46440", out$strain)
out$strain <- gsub("PA14_45800..PA14_46240.Δ40.707.b", "∆PA14_45800..PA14_46240", out$strain)
out$strain <- gsub("morA.E1153K.tRNAThr.tufB.G..A", "morA E1153K tRNAThr/tufB G->A", out$strain)
out$strain <- gsub("morA.N1124K", "morA N1124K", out$strain)
out$strain <- gsub("morA.K1123E", "morA K1123E", out$strain)

out$genotype <- "Ancestor"
out$genotype[grep("lasR", out$strain)] <- "lasR"
out$genotype[grep("PA14_", out$strain)] <- "lasR"
out$genotype[grep("morA", out$strain)] <- "morA"

######################

### k ###
out %>% 
  group_by(genotype, media) %>%
  summarise(mean = mean(k, na.rm = TRUE),
            sd = sd(k, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se) %>%
  ggplot(mapping = aes(x = media, y = mean, color = genotype)) +
  geom_point(position=position_dodge(.9), size =3) +
  theme_minimal() + theme(legend.position = "none") + ylab("Carrying Capacity (k)") +
  scale_color_manual(values = c("Ancestor" = "dark gray", "lasR" = "#4393C3", "morA" = "#339966")) + 
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.5, size =0.75,
                position=position_dodge(.9)) 
ggsave(filename= "k.png", dpi=300, dev='png', width = 3, height = 2)

out %>% 
  group_by(genotype, media) %>%
  ggplot(mapping = aes(x = media, y = k, color = genotype)) +
  geom_point(position=position_dodge(.9), size =3) +
  theme_minimal() + theme(legend.position = "none") + ylab("Carrying Capacity (k)") +
  scale_color_manual(values = c("Ancestor" = "dark gray", "lasR" = "#4393C3", "morA" = "#339966")) 

### ANOVA
k <- lm(k ~ media+genotype+media:genotype, data = out)
k.aov <- aov(k)
summary(k.aov)
### Tukey
letters.df <- data.frame(multcompLetters(TukeyHSD(k.aov)$'media:genotype'[,4])$Letters)

### auc_e ###
out %>% 
  group_by(genotype, media) %>%
  summarise(mean = mean(auc_e, na.rm = TRUE),
            sd = sd(auc_e, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se) %>%
  mutate(media, fct_relevel(media, "all", "glucose", "amino acids", "lactate")) %>%
  ggplot(mapping = aes(x = media, y = mean, color = genotype)) +
  geom_point(position=position_dodge(.9), size=3) +
  theme_minimal() + theme(legend.position = "none") + ylab("Area Under Curve") +
  scale_color_manual(values = c("Ancestor" = "dark gray", "lasR" = "#4393C3", "morA" = "#339966")) + 
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.5, size=0.75,
                position=position_dodge(.9)) 
ggsave(filename= "auc.png", dpi=300, dev='png', width = 3, height = 2)

auc <- lm(auc ~ media+genotype+media:genotype, data = out)
auc.aov <- aov(auc)
letters.df <- data.frame(multcompLetters(TukeyHSD(auc.aov)$'media:genotype'[,4])$Letters)
colnames(letters.df)[1] <- "Letter" 
letters.df$Category <- rownames(letters.df) 
letters.df <- letters.df[order(letters.df$Category) , ]

### r ###
out %>% 
  group_by(genotype, media) %>%
  summarise(mean = mean(r, na.rm = TRUE),
            sd = sd(r, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se) %>%
  mutate(media, fct_relevel(media, "all", "glucose", "amino acids", "lactate")) %>%
  ggplot(mapping = aes(x = media, y = mean, color = genotype)) +
  geom_point(position=position_dodge(.9), size=3) +
  theme_minimal() + theme(legend.position = "none") + ylab("Max Growth Rate (r)") +
  scale_color_manual(values = c("Ancestor" = "dark gray", "lasR" = "#4393C3", "morA" = "#339966")) + 
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.5, size=0.75,
                position=position_dodge(.9)) 
ggsave(filename= "legend.png", dpi=300, dev='png', width = 3, height = 2)

r <- lm(r ~ media+genotype+media:genotype, data = out)
r.aov <- aov(r)
letters.df <- data.frame(multcompLetters(TukeyHSD(r.aov)$'media:genotype'[,4])$Letters)
colnames(letters.df)[1] <- "Letter" 
letters.df$Category <- rownames(letters.df) 
letters.df <- letters.df[order(letters.df$Category) , ]
