setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(illuminaio)
library(gplots)
library(seqinr)
library(plyr)
library(grid)


#Downloaded from Table Browser
# Group: Genes and Gene Prediction
# Track: NCBI RefSeq
# Output Format: Selected Field from primary and related tables
hgTables <- read.delim("Data/hgTables Summary.txt")


genelist <- readBGX("Data/HumanHT-12_V3_0_R3_11283641_A.bgx")
genelist <- data.frame(genelist$probes)
genelist <- data.frame(Gene.Symbol=genelist$ILMN_Gene,
                       Probes=genelist$Probe_Id,
                       Refseq=genelist$RefSeq_ID,
                       Unigene=genelist$Unigene_ID,
                       Chromosome=genelist$Chromosome,
                       Coordinates=genelist$Probe_Coordinates)
genelist <- genelist[-grep("NR",genelist$Refseq),] # 654 genes
genelist <- genelist[-grep("XR",genelist$Refseq),] # 501 genes
genelist$Coordinates[genelist$Coordinates==""] <- NA
genelist <- genelist[complete.cases(genelist),] # 42084 probes w/ coordinates
rownames(genelist) <- genelist$Probes

hgTables <- hgTables[(na.omit(match(hgTables$name2,genelist$Gene.Symbol))),]


breast <- read.table("Data/Outputs/Breast P1.txt")
breast$Tissue <- "Breast"
cerebellum <- read.table("Data/Outputs/Cerebellum P1.txt")
cerebellum$Tissue <- "Cerebellum"
frontal <- read.table("Data/Outputs/Frontal P1.txt")
frontal$Tissue <- "Frontal Cortex"

df2 <- rbind(breast, cerebellum, frontal)
df2 <- merge(df2, genelist[,c(1,2)], by = "Probes")
df2$Class <- revalue(df2$Class, c("NV" = "Non-Variable"))

tissues <- c("Breast", "Cerebellum", "Frontal Cortex")
genesets <- c("Hypervariable", "Hypovariable", 'Non-Variable')


pdf("Figures/Additional File 1.pdf")
textplot("Additional File 1. Structural analysis of genes as a\nfunction of EV
", halign="center", valign='top')

######################
# Largest Transcript #
######################
tx.size <- unique(data.frame(Gene.Symbol=hgTables$name2,
                             tx.size=(hgTables$txEnd-hgTables$txStart)/1000))

df <- merge(df2, tx.size, by="Gene.Symbol")
df <- df[order(df$Gene.Symbol, -df$tx.size),]
df <- df[!duplicated(df$Gene.Symbol),]

lm_function <- function(df) {
  m <- lm(EV~tx.size, df)
  eq <- data.frame("Intercept" = coef(m)[[1]],
                   "Slope" = coef(m)[[2]],
                   "P-value" = summary(m)$coefficients[2,4])
}

eq <- ddply(df,.(Tissue, Class),lm_function)
eq_table <- grob(textplot(eq))
title("Largest Transcript Size Linear Regression Analysis")
eq_table

for (i in tissues) {
  for (j in genesets){
    plot_df <- subset(df, df$Tissue == i & df$Class == j)
    m <- lm(EV~tx.size, plot_df)
    
    plot_title <- paste(i,j)
    plot_text <- paste('Intercept: ',format(coef(m)[[1]], digits=3),
                       '   Slope: ', format(coef(m)[[2]], digits=3),
                       '   R2: ', format(summary(m)$adj.r.squared, digits=3),
                       '   Correlation: ', format(cor(plot_df$tx.size, plot_df$EV, method='kendall'), digits=3))
    plot(plot_df$tx.size, plot_df$EV, lwd=1, pch=16,
         xlab="Transcript Size\n", ylab='EV', main=plot_title, sub = plot_text)
    abline(a=coef(m)[[1]], b=coef(m)[[2]],col='red', lwd=2)
  }
}



#######################
# Smallest Transcript #
#######################
df <- merge(df2, tx.size, by="Gene.Symbol")
df <- df[order(df$Gene.Symbol, df$tx.size),]
df <- df[!duplicated(df$Gene.Symbol),]


lm_function <- function(df) {
  m <- lm(EV~tx.size, df)
  eq <- data.frame("Intercept" = coef(m)[[1]],
                   "Slope" = coef(m)[[2]],
                   "P-value" = summary(m)$coefficients[2,4])
}

eq <- ddply(df,.(Tissue, Class),lm_function)
eq_table <- grob(textplot(eq))
title("Smallest Transcript Size Linear Regression Analysis")
eq_table

for (i in tissues) {
  for (j in genesets){
    plot_df <- subset(df, df$Tissue == i & df$Class == j)
    m <- lm(EV~tx.size, plot_df)
    
    plot_title <- paste(i,j)
    plot_text <- paste('Intercept: ',format(coef(m)[[1]], digits=3),
                       '   Slope: ', format(coef(m)[[2]], digits=3),
                       '   R2: ', format(summary(m)$adj.r.squared, digits=3),
                       '   Correlation: ', format(cor(plot_df$tx.size, plot_df$EV, method='kendall'), digits=3))
    plot(plot_df$tx.size, plot_df$EV, lwd=1, pch=16,
         xlab="Transcript Size\n", ylab='EV', main=plot_title, sub = plot_text)
    abline(a=coef(m)[[1]], b=coef(m)[[2]],col='red', lwd=2)
  }
}

##############
# Exon Count #
##############
exon.count <- aggregate(exonCount~name2,data=hgTables,FUN=median)
colnames(exon.count) <- c("Gene.Symbol", "exon.count")
df <- merge(df2, exon.count, by="Gene.Symbol")


lm_function <- function(df) {
  m <- lm(EV~exon.count, df)
  eq <- data.frame("Intercept" = coef(m)[[1]],
                   "Slope" = coef(m)[[2]],
                   "P-value" = summary(m)$coefficients[2,4])
}

eq <- ddply(df,.(Tissue, Class),lm_function)
eq_table <- grob(textplot(eq))
title("Exon Count Linear Regression Analysis")
eq_table

for (i in tissues) {
  for (j in genesets){
    plot_df <- subset(df, df$Tissue == i & df$Class == j)
    m <- lm(EV~exon.count, plot_df)
    
    plot_title <- paste(i,j)
    plot_text <- paste('Intercept: ',format(coef(m)[[1]], digits=3),
                       '   Slope: ', format(coef(m)[[2]], digits=3),
                       '   R2: ', format(summary(m)$adj.r.squared, digits=3),
                       '   Correlation: ', format(cor(plot_df$exon.count, plot_df$EV, method='kendall'), digits=3))
    plot(plot_df$exon.count, plot_df$EV, lwd=1, pch=16,
         xlab="Transcript Size\n", ylab='EV', main=plot_title, sub = plot_text)
    abline(a=coef(m)[[1]], b=coef(m)[[2]],col='red', lwd=2)
  }
}



# Sequence Length
s.length <- read.fasta("Data/hgTable Sequences.txt")
s.length1 <- getLength(s.length)
s.length2 <- sapply(strsplit(getName(s.length),split="Seq_", fixed=TRUE), function(x) (x[2]))
s.length <- data.frame(Refseq=s.length2, Seq.len=s.length1)
s.length <- merge(genelist[,1:3,drop=F],s.length, by="Refseq")

df <- merge(df2, s.length, by="Gene.Symbol")

lm_function <- function(df) {
  m <- lm(EV~Seq.len, df)
  eq <- data.frame("Intercept" = coef(m)[[1]],
                   "Slope" = coef(m)[[2]],
                   "P-value" = summary(m)$coefficients[2,4])
}

eq <- ddply(df,.(Tissue, Class),lm_function)
eq_table <- grob(textplot(eq))
title("Sequence Length Linear Regression Analysis")
eq_table

for (i in tissues) {
  for (j in genesets){
    plot_df <- subset(df, df$Tissue == i & df$Class == j)
    m <- lm(EV~Seq.len, plot_df)
    
    plot_title <- paste(i,j)
    plot_text <- paste('Intercept: ',format(coef(m)[[1]], digits=3),
                       '   Slope: ', format(coef(m)[[2]], digits=3),
                       '   R2: ', format(summary(m)$adj.r.squared, digits=3),
                       '   Correlation: ', format(cor(plot_df$Seq.len, plot_df$EV, method='kendall'), digits=3))
    plot(plot_df$Seq.len, plot_df$EV, lwd=1, pch=16,
         xlab="Transcript Size\n", ylab='EV', main=plot_title, sub = plot_text)
    abline(a=coef(m)[[1]], b=coef(m)[[2]],col='red', lwd=2)
  }
}

#########################
# Number of Transcripts #
#########################

n.tx <- ddply(hgTables,.(name2),nrow)
colnames(n.tx) <- c("Gene.Symbol", "n.tx")
df <- merge(df2, n.tx, by="Gene.Symbol")


lm_function <- function(df) {
  m <- lm(EV~n.tx, df)
  eq <- data.frame("Intercept" = coef(m)[[1]],
                   "Slope" = coef(m)[[2]],
                   "P-value" = summary(m)$coefficients[2,4])
}

eq <- ddply(df,.(Tissue, Class),lm_function)
eq_table <- grob(textplot(eq))
title("Number of Transcripts Linear Regression Analysis")
eq_table

for (i in tissues) {
  for (j in genesets){
    plot_df <- subset(df, df$Tissue == i & df$Class == j)
    m <- lm(EV~n.tx, plot_df)
    
    plot_title <- paste(i,j)
    plot_text <- paste('Intercept: ',format(coef(m)[[1]], digits=3),
                       '   Slope: ', format(coef(m)[[2]], digits=3),
                       '   R2: ', format(summary(m)$adj.r.squared, digits=3),
                       '   Correlation: ', format(cor(plot_df$n.tx, plot_df$EV, method='kendall'), digits=3))
    plot(plot_df$n.tx, plot_df$EV, lwd=1, pch=16,
         xlab="Transcript Size\n", ylab='EV', main=plot_title, sub = plot_text)
    abline(a=coef(m)[[1]], b=coef(m)[[2]],col='red', lwd=2)
  }
}



dev.off()