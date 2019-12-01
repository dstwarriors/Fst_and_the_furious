#packages
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Install and load libraries
packages(tidyr)
packages(dplyr)
packages(foreach)
packages(BSDA)
packages(ggplot2)
packages(tidyverse)
packages(openxlsx)
library(tidyr)
library(dplyr)
library(foreach)
library(BSDA)
library(ggplot2)
library(tidyverse)
library(openxlsx)

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
topSNPs <- read.xlsx("41588_2018_147_MOESM3_ESM.xlsx", sheet=2, startRow=1, colNames=TRUE)
TenkSNPs <- read.delim("GWAS_EA.to10K.txt", header=TRUE, stringsAsFactors=FALSE)

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
colnames(topSNPs) <- c("SNP", "CHR", "BP", "ALL1", "ALL2", "Frq1_TOP", "Effect_Size_Top", "SE_TOP", "PVal_Top", "HeteroI", "HeteroP", "N_TOP")
colnames(TenkSNPs) <- c("SNP", "CHR", "BP", "ALLE1", "ALLE2", "EU_FRQ", "BETA_TEN", "SE_TEN", "PVAL_TEN")
#Pick columns

WF_GWAS_SNPs <- WF_GWAS_SNPs[,names(WF_GWAS_SNPs) %in% c("MARKER_ADDED", "ALLELE_1", "ALLELE_2", "Freq1", "Effect", "P_VALUE")]
AfrFrq <- AfrFrq[,names(AfrFrq) %in% c("CHR", "BP", "A1", "A1_Frq_Afr", "A2", "A2_Frq_Afr")]
EurFrq <- EurFrq[,names(EurFrq) %in% c("CHR", "BP", "A1", "A1_Frq_Eur", "A2", "A2_Frq_Eur")]
topSNPs <- topSNPs[,names(topSNPs) %in% c("CHR", "BP", "ALL1", "ALL2", "Frq1_TOP", "Effect_Size_Top", "PVal_Top")]
TenkSNPs <- TenkSNPs[,names(TenkSNPs) %in% c("CHR", "BP", "ALLE1", "ALLE2", "BETA_TEN", "PVAL_TEN")]

# Regex editing

WF_GWAS_SNPs$BP <- gsub("\\d+:", "", WF_GWAS_SNPs$MARKER_ADDED)
WF_GWAS_SNPs$CHR <- gsub(":\\d+", "", WF_GWAS_SNPs$MARKER_ADDED)
WF_GWAS_SNPs$ALLELE_1 <- str_to_upper(WF_GWAS_SNPs$ALLELE_1, locale = "en")
WF_GWAS_SNPs$ALLELE_2 <- str_to_upper(WF_GWAS_SNPs$ALLELE_2, locale = "en")

# Convert data type

AfrFrq$CHR <- as.character(AfrFrq$CHR)
EurFrq$CHR <- as.character(EurFrq$CHR)
BF_GWAS_SNPs$CHR <- as.character(BF_GWAS_SNPs$CHR)
topSNPs$CHR <- as.character(topSNPs$CHR)
TenkSNPs$CHR <- as.character(TenkSNPs$CHR)
AfrFrq$BP <- as.character(AfrFrq$BP)
EurFrq$BP <- as.character(EurFrq$BP)
BF_GWAS_SNPs$BP <- as.character(BF_GWAS_SNPs$BP)
topSNPs$BP <- as.character(topSNPs$BP)
TenkSNPs$BP <- as.character(TenkSNPs$BP)

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

# Join 

FinalTop <- inner_join(Final, topSNPs, by=c("CHR", "BP"))
FinalTen <- inner_join(Final, TenkSNPs, by=c("CHR", "BP"))

# Make sure allele ordering is proper

TF <- FinalTop$A1!=FinalTop$ALL1
Ha <- FinalTop$ALL1[TF]
Hb <- FinalTop$ALL2[TF]
Hc <- FinalTop$Effect_Size_top[TF]
FinalTop$ALL1[TF] <- Hb
FinalTop$ALL2[TF] <- Ha
FinalTop$Effect_Size_Top[TF] <- -1*Hc

TF <- FinalTen$A1!=FinalTen$ALLE1
Ja <- FinalTen$ALLE1[TF]
Jb <- FinalTen$ALLE2[TF]
Jc <- FinalTen$BETA_TEN[TF]
FinalTen$ALLE1[TF] <- Jb
FinalTen$ALLE2[TF] <- Ja
FinalTen$BETA_TEN[TF] <- -1*Jc

# Test if allele ordering is proper

C <- FinalTop$A1==FinalTop$ALL1
print(nrow(FinalTop)-length(C[C]))
D <- FinalTen$A1==FinalTen$ALLE1
print(nrow(FinalTen)-length(D[D]))

# Create diff

Final$Diff <- Final$A1_Frq_Afr - Final$A1_Frq_Eur
FinalTop$Diff <- FinalTop$A1_Frq_Afr-FinalTop$A1_Frq_Eur
FinalTen$Diff <- FinalTen$A1_Frq_Afr-FinalTen$A1_Frq_Eur

# Remove NAs

Final <- Final[complete.cases(Final$Effect), ]
Final <- Final[complete.cases(Final$Diff), ]
FinalTop <- FinalTop[complete.cases(FinalTop$Diff), ]

# Print number of rows

nrow(Final)

# Print results

num1 <- sum(Final$A1_Frq_Afr*Final$Effect)
num2 <- sum(Final$A1_Frq_Eur*Final$Effect)
num3 <- sum(Final$Diff*Final$Beta)
print(num1)
print(num2)

# Test against null

NullAllWF<-1:10000
NullAllBF<-1:10000
TempDF <- Final[,names(Final) %in% c("Effect", "Beta", "Diff")]
for(i in 1:10000){
	F<-TempDF
	Thing <- sample(1:nrow(TempDF),size=floor(nrow(TempDF)/2))
	F$Effect[Thing] <- -1*TempDF$Effect[Thing]
	F$Beta[Thing] <- -1*TempDF$Beta[Thing]
	NullAllWF[i] <- sum(F$Effect*F$Diff)^2
	NullAllBF[i] <- sum(F$Beta*F$Diff)^2
}
NullAllWF <- as.data.frame(NullAllWF)
NullAllBF <- as.data.frame(NullAllBF)
sigTestAllWF <-length(NullAllWF[NullAllWF>(num1-num2)^2])/nrow(NullAllWF)
sigTestAllBF <-length(NullAllBF[NullAllBF>(num3)^2])/nrow(NullAllBF)
colnames(NullAllBF) <- c("N")
colnames(NullAllWF) <- c("N")
NullAllWFPlot <- ggplot(data=NullAllWF, aes(x=N))+
	labs(x="PGS Difference Squared", title="EduYear ALL SNPs (Lee et. al 2018) - WF effect sizes")+
	geom_histogram(bins=50, fill="#377EB8", alpha=0.5)+
	annotate(geom="text", x=1.3, y=200, label=paste0("P=", round(sigTestAllWF, digits=3)))+
	geom_vline(xintercept = (num1-num2)^2)

NullAllBFPlot <- ggplot(data=NullAllBF, aes(x=N))+
	labs(x="PGS Difference Squared", title="EduYear ALL SNPs (Lee et. al 2018) - BF effect sizes")+
	geom_histogram(bins=50, fill="#377EB8", alpha=0.5)+
	annotate(geom="text", x=2, y=300, label=paste0("P=", round(sigTestAllBF, digits=3)))+
	geom_vline(xintercept = (num3)^2)
multiplot(NullAllBFPlot, NullAllWFPlot, cols=2)

readline(prompt="Press [enter] to continue")



# Only genome-wide significant alleles

df_sig5e8 <- Final[Final$Pval < 5e-8, ]

nrow(df_sig5e8)

num5 <- sum(df_sig5e8$A1_Frq_Afr*df_sig5e8$Effect)
num6 <- sum(df_sig5e8$A1_Frq_Eur*df_sig5e8$Effect)
num7 <- sum(df_sig5e8$A1_Frq_Afr*df_sig5e8$Beta)
num8 <- sum(df_sig5e8$A1_Frq_Eur*df_sig5e8$Beta)
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
NullBF <- as.data.frame(NullBF)
NullWF <- as.data.frame(NullWF)
sigTestWF <- length(NullWF[NullWF>((num5-num6)^2)])/nrow(NullWF)
sigTestBF <- length(NullBF[NullBF>((num7-num8)^2)])/nrow(NullBF)
colnames(NullBF) <- c("N")
colnames(NullWF) <- c("N")
NullWFPlot <- ggplot(data=NullWF, aes(x=N))+
	labs(x="PGS Difference Squared", title="EduYear GWS SNps (Lee et. al 2018) - WF effect sizes")+
	geom_histogram(bins=50, fill="#377EB8", alpha=0.5)+
	annotate(geom="text", x=0.015, y=500, label=paste0("P=", round(sigTestWF, digits=3)))+
	geom_vline(xintercept = (num5-num6)^2)


NullBFPlot <- ggplot(data=NullBF, aes(x=N))+
	labs(x="PGS Difference Squared", title="EduYear GWS SNPs (Lee et. al 2018) - BF effect sizes")+
	geom_histogram(bins=50, fill="#377EB8", alpha=0.5)+
	annotate(geom="text", x=0.025, y=500, label=paste0("P=", round(sigTestBF, digits=3)))+
	geom_vline(xintercept = (num7-num8)^2)

multiplot(NullBFPlot, NullWFPlot, cols=2)

readline(prompt="Press [enter] to continue")

# Using "top SNPs"

nrow(FinalTop)

num9 <- sum(FinalTop$A1_Frq_Afr*FinalTop$Effect_Size_Top)
num10 <- sum(FinalTop$A1_Frq_Eur*FinalTop$Effect)
num11 <- sum(FinalTop$A1_Frq_Afr*FinalTop$Beta)
num12 <- sum(FinalTop$A1_Frq_Eur*FinalTop$Beta)
print(num9)
print(num10)
print(num11)
print(num12)
NullWFTop<-1:10000
NullBFTop<-1:10000
TempDFTop <- FinalTop[,names(FinalTop) %in% c("Effect", "Beta", "Diff")]
for(i in 1:10000){
	F3<-TempDFTop
	ThingTop <- sample(1:nrow(TempDFTop),size=floor(nrow(TempDFTop)/2))
	F3$Effect[ThingTop] <- -1*TempDFTop$Effect[ThingTop]
	F3$Beta[ThingTop] <- -1*TempDFTop$Beta[ThingTop]
	NullWFTop[i] <- sum(F3$Effect*F3$Diff)^2
	NullBFTop[i] <- sum(F3$Beta*F3$Diff)^2
}
sigTestWFTop <- length(NullWFTop[NullWFTop>((num9-num10)^2)])/length(NullWFTop)
sigTestBFTop <- length(NullBFTop[NullBFTop>((num11-num12)^2)])/length(NullBFTop)
NullBFTop <- as.data.frame(NullBFTop)
NullWFTop <- as.data.frame(NullWFTop)
colnames(NullBFTop) <- c("N")
colnames(NullWFTop) <- c("N")
NullWFTopPlot <- ggplot(data=NullWFTop, aes(x=N))+
	labs(x="PGS Difference Squared", title="EduYear Top SNPs (Lee et. al 2018) - WF effect sizes")+
	geom_histogram(bins=50, fill="#377EB8", alpha=0.5)+
	annotate(geom="text", x=0.015, y=500, label=paste0("P=", round(sigTestWFTop, digits=3)))+
	geom_vline(xintercept = (num9-num10)^2)

NullBFTopPlot <- ggplot(data=NullBFTop, aes(x=N))+
	labs(x="PGS Difference Squared", title="EduYear Top SNPs (Lee et. al 2018) - BF effect sizes")+
	geom_histogram(bins=50, fill="#377EB8", alpha=0.5)+
	annotate(geom="text", x=0.025, y=500, label=paste0("P=", round(sigTestBFTop, digits=3)))+
	geom_vline(xintercept = (num11-num12)^2)
multiplot(NullBFTopPlot, NullWFTopPlot, cols=2)


Matr <- 1:10
for(i in 1:10){
	Tem <- Final[Final$Pval < 5*10^(-i),]
	Matr[i] <- nrow(Tem[(Tem$Effect>0 & Tem$Beta>0)|(Tem$Effect<0 & Tem$Beta<0),])/nrow(Tem)
}
print(Matr)
cor.test(rank(Matr), Matr)
cor.test(5*10^(-rank(Matr)), Matr)
