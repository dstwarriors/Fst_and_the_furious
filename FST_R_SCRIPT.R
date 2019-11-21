#################################
# DST WARRIORS
# FST and the furious code
#
# Written by Kevin Bird
# Tidied by freerecall
#
#################################

# Install and load libraries
install.packages("tidyr")
install.packages("dplyr")
install.packages("foreach")
install.packages("BSDA")
install.packages("ggplot2")
library(tidyr)
library(dplyr)
library(foreach)
library(BSDA)
library(ggplot2)

# Set working directory
setwd("D:/DST Warriors/Fst_Files")

# Set seed
set.seed(5)

# Read in dfs
all_phase3_MAF <- read.delim("all_phase3_MAF.txt", stringsAsFactors=FALSE)
EUR_LDSC_all <- read.delim("EUR_LDSC_all.txt", stringsAsFactors=FALSE)
AFR_EUR_Genomewide_Filtered.weir <- read.delim("AFR_EUR_Genomewide_Filtered.weir.fst")
GWAS_SNPs <- read.delim("GWAS_EA_SNPs_clumped5e6.txt", header=TRUE, stringsAsFactors=FALSE)

# Standardize col names
colnames(AFR_EUR_Genomewide_Filtered.weir) <- c("CHR", "BP", "FST")
colnames(GWAS_SNPs) <- "SNP"
head(AFR_EUR_Genomewide_Filtered.weir)
head(GWAS_SNPs)

# Combine dataframes
EUR_LDSC_all<-left_join(EUR_LDSC_all,all_phase3_MAF,by="SNP")
EUR_LDSC_all<-left_join(EUR_LDSC_all,AFR_EUR_Genomewide_Filtered.weir,by=c("CHR","BP"))

# Grab MAF
MAF <- all_phase3_MAF$MAF

# Calculate MAF Quintile and L2 quintile
ALLSNPs <- within(EUR_LDSC_all, MAF_quintile <- as.integer(cut(MAF, quantile(MAF, probs=0:20/20), include.lowest=TRUE)))
ALLSNPs <- within(ALLSNPs, L2_quintile <- as.integer(cut(L2, quantile(L2, probs=0:20/20), include.lowest=TRUE)))

# Combine df by SNP
GWAS_SNPs<-left_join(GWAS_SNPs,ALLSNPs,by="SNP")

# Filt out non-GWAS SNPS
NOGWASSNPs <-ALLSNPs %>% filter(!(SNP %in% GWAS_SNPs$SNP))

# Bootstrapping: THIS WILL TAKE AWHILE
NullFst<-as.data.frame(foreach(row=1:nrow(GWAS_SNPs), .combine='rbind') %do% 
                        sample_n(NOGWASSNPs %>% 
                        filter(MAF_quintile==GWAS_SNPs$MAF_quintile[row] & L2_quintile==GWAS_SNPs$L2_quintile[row]),10000,replace = T)$FST)


# Set negative and null values to 0
NullFst[NullFst < 0]<-0
NullFst[is.na(NullFst)]<-0
GWAS_SNPs$FST[GWAS_SNPs$FST<0]<-0

# Compute means
NullFstDist <- NullFst %>% 
  summarise_each(funs(mean)) %>% 
  gather() %>% 
  select(-key)

# Signiifiance testing
Sigtest<-z.test(GWAS_SNPs$FST,alternative="two.sided",mu=mean(NullFstDist$value),sigma.x = sd(GWAS_SNPs$FST))

# Plot bootstrapped means
ggplot(data=NullFstDist,aes(x=value))+
  labs(x="Mean Fst",title = "EduYear (Lee et al. 2018)")+
  geom_histogram(bins=50,fill="#377EB8",alpha=0.5)+
  annotate(geom="text",x=0.0655,y=500, label=paste0("P=",round(Sigtest$p.value,digits = 3)))+
  geom_vline(xintercept = mean(GWAS_SNPs$FST))
