#packages
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

# Install and load libraries
packages(tidyr)
packages(dplyr)
packages(foreach)
packages(BSDA)
packages(ggplot2)
packages(tidyverse)
library(tidyr)
library(dplyr)
library(foreach)
library(BSDA)
library(ggplot2)
library(tidyverse)

# Set working directory
setwd("C:/Users/<>/Documents/Fst/1KG")

# Set seed
set.seed(5)

# Read in dfs
WF_GWAS_SNPs <- read.delim("WF_STR_UKB_WLS_clumped_meta1.TBL", header=TRUE, stringsAsFactors=FALSE)
AfrFrq <- read.delim("Lee_etal_WithinFam_AfrFreq.frq", header=TRUE, stringsAsFactors=FALSE)
EurFrq <- read.delim("Lee_etal_WithinFam_EurFreq.frq", header=TRUE, stringsAsFactors=FALSE)
BF_GWAS_SNPs <- read.delim("GWAS_EA_excl23andMe.txt", header=TRUE, stringsAsFactors=FALSE)
LDSC <- read.delim("EUR_LDSC_all.txt", header=TRUE, stringsAsFactors=FALSE)

#Pick columns
BF_GWAS_SNPs <- BF_GWAS_SNPs[,names(BF_GWAS_SNPs) %in% c("CHR", "POS", "A1", "A2", "Beta", "Pval")]

#Print

head(WF_GWAS_SNPs)
head(AfrFrq)
head(EurFrq)

# Column names!!
# Standardize col names

colnames(WF_GWAS_SNPs) <- c("MARKER_ADDED", "ALLELE_1", "ALLELE_2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "STD_ERR", "P_VALUE", "DIRECTION", "SAMPLE_SIZE")
colnames(AfrFrq) <- c("CHR", "BP", "N_ALLELES", "N_CHR_Afr", "A1", "A1_Frq_Afr", "A2", "A2_Frq_Afr")
colnames(EurFrq) <- c("CHR", "BP", "N_ALLELES", "N_CHR_Eur", "A1", "A1_Frq_Eur", "A2", "A2_Frq_Eur")
colnames(BF_GWAS_SNPs) <- c("CHR", "BP", "AL1", "AL2", "Beta", "Pval")

#Pick columns

WF_GWAS_SNPs <- WF_GWAS_SNPs[,names(WF_GWAS_SNPs) %in% c("MARKER_ADDED", "ALLELE_1", "ALLELE_2", "Freq1", "Effect", "P_VALUE")]
AfrFrq <- AfrFrq[,names(AfrFrq) %in% c("CHR", "BP", "A1", "A1_Frq_Afr", "A2", "A2_Frq_Afr")]
EurFrq <- EurFrq[,names(EurFrq) %in% c("CHR", "BP", "A1", "A1_Frq_Eur", "A2", "A2_Frq_Eur")]

# Regex editing

WF_GWAS_SNPs$BP <- gsub("\\d+:", "", WF_GWAS_SNPs$MARKER_ADDED)
WF_GWAS_SNPs$CHR <- gsub(":\\d+", "", WF_GWAS_SNPs$MARKER_ADDED)
WF_GWAS_SNPs$ALLELE_1 <- str_to_upper(WF_GWAS_SNPs$ALLELE_1, locale = "en")
WF_GWAS_SNPs$ALLELE_2 <- str_to_upper(WF_GWAS_SNPs$ALLELE_2, locale = "en")

# Convert data type

AfrFrq$CHR <- as.character(AfrFrq$CHR)
EurFrq$CHR <- as.character(EurFrq$CHR)
BF_GWAS_SNPs$CHR <- as.character(BF_GWAS_SNPs$CHR)
AfrFrq$BP <- as.character(AfrFrq$BP)
EurFrq$BP <- as.character(EurFrq$BP)
BF_GWAS_SNPs$BP <- as.character(BF_GWAS_SNPs$BP)

# Join

Frqs <- inner_join(AfrFrq, EurFrq, by=c("CHR", "BP", "A1", "A2"))
Im <- inner_join(Frqs, WF_GWAS_SNPs, by=c("CHR", "BP"))
Final <- inner_join(Im, BF_GWAS_SNPs, by=c("CHR", "BP"))

# Make sure the allele ordering is correct

TF <- Final$A1!=Final$ALLELE_1
K <- Final$A1[TF]
L <- Final$A2[TF]
Final$A1[TF] <- L
Final$A2[TF] <- K
Final$Freq1[TF] <- 1-Final$Freq1[TF]
K <- Final$A1_Frq_Afr[TF]
L <- Final$A2_Frq_Afr[TF]
Final$A1_Frq_Afr[TF] <- L
Final$A2_Frq_Afr[TF] <- K
K <- Final$A1_Frq_Eur[TF]
L <- Final$A2_Frq_Eur[TF]
Final$A1_Frq_Eur[TF] <- L
Final$A2_Frq_Eur[TF] <- K

FT <- Final$A1!=Final$AL1
K <- Final$AL1[FT]
L <- Final$AL2[FT]
M <- Final$Beta[FT]
Final$AL1[FT] <- L
Final$AL2[FT] <- K
Final$Beta[FT] <- -1*M

#Test if primary alleles were properly matched

A<-Final$A1==Final$ALLELE_1
print(nrow(Final)-length(A[A]))
B<-Final$A1==Final$AL1
print(nrow(Final)-length(B[B]))

# Create diff

Final$Diff <- Final$A1_Frq_Afr - Final$A1_Frq_Eur

# Remove NAs

Final <- Final[complete.cases(Final$Effect), ]
Final <- Final[complete.cases(Final$Diff), ]


# Print number of rows

nrow(Final)

# Print results

num1 <- sum(Final$A1_Frq_Afr*Final$Effect)
num2 <- sum(Final$A1_Frq_Eur*Final$Effect)
num3 <- sum(Final$Diff*Final$Beta)
print(num1)
print(num2)

# Test against null

numm <- 1000
K<-1:numm
L<-1:numm
TempDF <- Final[,names(Final) %in% c("Effect", "Beta", "Diff")]
for(i in 1:numm){
	F<-TempDF
	Thing <- sample(1:nrow(TempDF),size=floor(nrow(TempDF)/2))
	F$Effect[Thing] <- -1*TempDF$Effect[Thing]
	F$Beta[Thing] <- -1*TempDF$Beta[Thing]
	K[i] <- sum(F$Effect*(F$A1_Frq_Afr-F$A1_Frq_Eur))^2
	L[i] <- sum(F$Beta*(F$A1_Frq_Afr-F$A1_Frq_Eur))^2
}
hist(K)
hist(L)
print(length(K[K>(num1-num2)^2])/length(K))
print(length(L[L>(num3)^2])/length(L))

# Only genome-wide significant alleles

df_sig5e8 <- Final[Final$Pval < 5e-8, ]

nrow(df_sig5e8)

num5 <- sum(df_sig5e8$A1_Frq_Afr*FS$Effect)
num6 <- sum(df_sig5e8$A1_Frq_Eur*FS$Effect)
num7 <- sum(df_sig5e8$A1_Frq_Afr*FS$Beta)
num8 <- sum(df_sig5e8$A1_Frq_Eur*FS$Beta)
print(num5)
print(num6)
print(num7)
print(num8)

NullWF<-1:10000
NullBF<-1:10000
TempDF <- df_sig5e8[,names(df_sig5e8) %in% c("Effect", "Beta", "Diff")]
for(i in 1:10000){
	F2<-TempDF
	Thing <- sample(1:nrow(TempDF),size=floor(nrow(TempDF)/2))
	F2$Effect[Thing] <- -1*TempDF$Effect[Thing]
	F2$Beta[Thing] <- -1*TempDF$Beta[Thing]
	NullWF[i] <- sum(F2$Effect*F2$Diff)^2
	NullBF[i] <- sum(F2$Beta*F2$Diff)^2
}
hist(NullWF)
hist(NullBF)
print(length(NullWF[NullWF>((num5-num6)^2)])/length(NullWF))
print(length(NullBF[NullBF>((num7-num8)^2)])/length(NullBF))

Matr <- 1:10
for(i in 1:10){
	Tem <- Final[Final$Pval < 5*10^(-i),]
	Matr[i] <- nrow(Tem[(Tem$Effect>0 & Tem$Beta>0)|(Tem$Effect<0 & Tem$Beta<0),])/nrow(Tem)
}
