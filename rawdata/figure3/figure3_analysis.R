### Analysis of motility and biofilm data for P. aerugionsa clones

library(readxl)
library(Rmisc)
library(tidyverse)
library(emmeans)
library(multcomp)
library(car)

setwd("/Users/mrs/Documents/pa14_nodrug/submit/rawdata/figure3")

### Swimming Motility ###

data <- read_excel("figure3_swimming.xlsx",sheet = 1)
data <- data[, 2:ncol(data)]
data <- pivot_longer(data = data, cols = 1:ncol(data), names_to = "genotype", values_to = "motility", values_drop_na = T)

sum <- summarySE(data,
                measurevar="motility",
                groupvars=c("genotype"))
pd = position_dodge(.2)

ggplot(sum, aes(x=genotype,y=motility)) +
  geom_errorbar(aes(ymin=motility-se,ymax=motility+se),width=.2, size=0.7, position=pd) +
  geom_point(shape=15, size=4, position=pd) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"), axis.text.x = element_text(angle = 90, hjust = 0.95)) 

data$genotype <- as.factor(data$genotype)
m <- lm(motility ~ genotype,data=data)
model <- aov(m)
summary(model)

#check assumptions
hist(residuals(model),
     col="darkgray") #histogram of residuals should appear normal
plot(fitted(model),
     residuals(model))
# homogeneity of variance
leveneTest(motility ~ genotype, data = data)
# normality
# Extract the residuals
aov_residuals <- residuals(object = model )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

# multiple comparisons
tuk <- glht(model, linfct = mcp(genotype = "Tukey"))
plot(tuk)
tuk.cld <- cld(tuk)
tuk.cld

#############################################################

### Biofilm ###

data <- read_excel("figure3_biofilm.xlsx",sheet = 1)
data <- data[, 2:ncol(data)]
data <- pivot_longer(data = data, cols = 1:ncol(data), names_to = "genotype", values_to = "biofilm", values_drop_na = T)

sum <- summarySE(data,
                 measurevar="biofilm",
                 groupvars=c("genotype"))
pd = position_dodge(.2)

ggplot(sum, aes(x=genotype,y=biofilm)) +
  geom_errorbar(aes(ymin=biofilm-se,ymax=biofilm+se),width=.2, size=0.7, position=pd) +
  geom_point(shape=15, size=4, position=pd) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"), axis.text.x = element_text(angle = 90, hjust = 0.95)) 

data$genotype <- as.factor(data$genotype)
m <- lm(biofilm ~ genotype,data=data)
model <- aov(m)
summary(model)

#check assumptions
hist(residuals(model),
     col="darkgray") #histogram of residuals should appear normal
plot(fitted(model),
     residuals(model))
# homogeneity of variance
leveneTest(biofilm ~ genotype, data = data)
# normality
# Extract the residuals
aov_residuals <- residuals(object = model )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

# multiple comparisons
tuk <- glht(model, linfct = mcp(genotype = "Tukey"))
plot(tuk)
tuk.cld <- cld(tuk)
tuk.cld

