setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(reshape2)
library(illuminaio)
library(limma)
library(dplyr)
library(ggplot2)
library(viridis)
library(plyr)
library(mixtools)



# Batch Effect Correction
meth.annotations <- read.table("Data/Methylation Annotations.txt",sep="\t",header=T)
colnames(meth.annotations) <- as.character(unlist(meth.annotations[2,]))
meth.annotations <- meth.annotations[c(1,2,10,11,12,13,14,15),-1]
meth.annotations <- as.data.frame(t(meth.annotations))
colnames(meth.annotations) <- c("Title", "Accession", "Tissue", "Sex", "Age", "PMI", "Tissuebank", "Batch")
samples <- rownames(meth.annotations)

meth.annotations <- apply(meth.annotations, 2, function(dat) {
  dat <- as.character(dat)
  dat <- gsub(".*: ","", dat)
})

meth.annotations <- as.data.frame(meth.annotations)
rownames(meth.annotations) <- samples
meth.annotations$Batch.Index <- as.factor(paste0(meth.annotations$Tissuebank, meth.annotations$Batch))

meth <- read.table("Data/Brain Methylation.txt", sep="\t")
rownames(meth) <- meth[,1]
IDs <- as.character(unlist(unname(meth[2,])))
IDs <- IDs[-1]
meth <- meth[-c(1,2),-1]
colnames(meth) <- IDs

indx <- sapply(meth, is.factor)
meth[indx] <- lapply(meth[indx], function(x) as.numeric(as.character(x)))

meth.annotations <- meth.annotations[order(match(meth.annotations$Accession,colnames(meth))),]
meth <- removeBatchEffect(meth, meth.annotations$Batch.Index)

# Match Gene & Methylation
genelist <- readBGX("Data/HumanHT-12_V3_0_R3_11283641_A.bgx")
genelist <- data.frame(genelist$probes)
genelist <- data.frame(Gene.Symbol=genelist$ILMN_Gene,
                       Probes=genelist$Probe_Id)
meth_probes <- read.csv("Data/Methylation Probes.csv")
meth_probes <- data.frame("ID"=meth_probes$ID, "Gene.Symbol"=meth_probes$Symbol)
meth_probes <- merge(meth_probes, genelist, by="Gene.Symbol")

# Match samples
exp_samples <- read.csv("Data/Brain Samples.csv",row.names=NULL, header=F)
exp_samples <- as.data.frame(t(exp_samples[1:2,]))
colnames(exp_samples) <- c("Sample", "Exp.Sample")

meth_samples <- meth.annotations[,1:2]
colnames(meth_samples) <- c("Sample", "Meth.Sample")
meth_samples$Sample <- gsub("-CpG", "", meth_samples$Sample)

sampledf <- merge(exp_samples, meth_samples, by="Sample")

df <- read.table("Data/Outputs/Final All 3.txt")
cere.sig <- subset(df, df$Tissue == "Cerebellum" & df$Sig == "Sig")
frontal.sig <- subset(df, df$Tissue == "Frontal Cortex" & df$Sig == "Sig")

cere <- read.table("Data/Outputs/Cerebellum Unimodal.txt")
cere <- cere[na.omit(match(cere.sig$Probes, rownames(cere))),
             na.omit(match(sampledf$Exp.Sample, colnames(cere)))]
cere$Probes <- rownames(cere)
cere <- merge(cere, meth_probes[,c(1,3)], by="Probes")
cere <- cere[ , -which(names(cere) %in% c("Probes"))]
cere <- melt(cere, id.vars = "Gene.Symbol")
colnames(cere) <- c("Gene.Symbol", "Sample", "Expression")
cere$Tissue <- "Cerebellum"
cere <- merge(cere, cere.sig[,c(2,4)], by="Gene.Symbol")
cere <- unique(cere)

frontal <- read.table("Data/Outputs/Frontal Unimodal.txt")
frontal <- frontal[na.omit(match(frontal.sig$Probes, rownames(frontal))),
             na.omit(match(sampledf$Exp.Sample, colnames(frontal)))]
frontal$Probes <- rownames(frontal)
frontal <- merge(frontal, meth_probes[,c(1,3)], by="Probes")
frontal <- frontal[ , -which(names(frontal) %in% c("Probes"))]
frontal <- melt(frontal, id.vars = "Gene.Symbol")
colnames(frontal) <- c("Gene.Symbol", "Sample", "Expression")
frontal$Tissue <- "Frontal Cortex"
frontal <- merge(frontal, frontal.sig[,c(2,4)], by="Gene.Symbol")
frontal <- unique(frontal)


# melt and add the methylation data
meth <- as.data.frame(meth)
mcere_samples <- sampledf[unique(na.omit(match(cere$Sample, sampledf$Exp.Sample))),]
mcere_probes <- meth_probes[na.omit(match(cere$Gene.Symbol, meth_probes$Gene.Symbol)),]
cere_meth <- meth[na.omit(match(unique(mcere_probes$ID), rownames(meth))),
                  na.omit(match(unique(mcere_samples$Meth.Sample), colnames(meth)))]
cere_meth$ID <- rownames(cere_meth)
cere_meth <- merge(cere_meth, meth_probes[,1:2], by="ID")
cere_meth <- cere_meth[ , -which(names(cere_meth) %in% c("ID"))]
cere_meth <- melt(cere_meth)
cere_meth <- unique(cere_meth)
colnames(cere_meth) <- c("Gene.Symbol", "Meth.Sample", "Methylation")
cere_meth <- merge(cere_meth, sampledf, by="Meth.Sample")
colnames(cere_meth)[5] <- "Sample"
cere_df <- merge(cere, cere_meth[,c(2,3,5)], by=c("Gene.Symbol", "Sample"))


mfrontal_samples <- sampledf[unique(na.omit(match(frontal$Sample, sampledf$Exp.Sample))),]
mfrontal_probes <- meth_probes[na.omit(match(frontal$Gene.Symbol, meth_probes$Gene.Symbol)),]
frontal_meth <- meth[na.omit(match(unique(mfrontal_probes$ID), rownames(meth))),
                  na.omit(match(unique(mfrontal_samples$Meth.Sample), colnames(meth)))]
frontal_meth$ID <- rownames(frontal_meth)
frontal_meth <- merge(frontal_meth, meth_probes[,1:2], by="ID")
frontal_meth <- frontal_meth[ , -which(names(frontal_meth) %in% c("ID"))]
frontal_meth <- melt(frontal_meth)
frontal_meth <- unique(frontal_meth)
colnames(frontal_meth) <- c("Gene.Symbol", "Meth.Sample", "Methylation")
frontal_meth <- merge(frontal_meth, sampledf, by="Meth.Sample")
colnames(frontal_meth)[5] <- "Sample"
frontal_df <- merge(frontal, frontal_meth[,c(2,3,5)], by=c("Gene.Symbol", "Sample"))

df <- rbind(frontal_df, cere_df)
df <- df[complete.cases(df),]


dat <- ddply(df, .(Gene.Symbol, Tissue), summarise, 
             Meth.Median = median(Methylation, na.rm = T),
             Exp.Median = median(Expression, na.rm = T))


cor.func <- function(x){
  x$Tissue <- as.character(x$Tissue)
  y1 <- subset(x, x$Tissue == "Cerebellum")
  x1.vals <- do.call(rbind,by(y1, y1$Gene.Symbol, function(x){
    x1_test <- cor.test(x$Expression, x$Methylation, method='kendall', exact = F)
    vals1 <- data.frame("Tissue" = "Cerebellum",
                        "pval" = x1_test$p.value,
                        "COR" = x1_test$estimate[[1]])
    vals1
  }))
  
  y2 <- subset(x, x$Tissue == "Frontal Cortex")
  x2.vals <- do.call(rbind,by(y2, y2$Gene.Symbol, function(x){
    x2_test <- cor.test(x$Expression, x$Methylation, method='kendall', exact = F)
    vals2 <- data.frame("Tissue" = "Frontal Cortex",
                        "pval" = x2_test$p.value,
                        "COR" = x2_test$estimate[[1]])
    vals2
  }))
  
  x1.vals$Gene.Symbol <- rownames(x1.vals)
  x2.vals$Gene.Symbol <- rownames(x2.vals)
  
  vals <- rbind(x1.vals, x2.vals)
  vals
}


paired.cor <- cor.func(df)

dat <- merge(dat, paired.cor, by=c("Gene.Symbol", "Tissue"), all=TRUE)
dat <- merge(df[,c(1,4,5)], dat, by=c("Gene.Symbol", "Tissue"))
dat <- unique(dat)



cere.dat <- subset(dat, dat$Tissue == "Cerebellum")
frontal.dat <- subset(dat, dat$Tissue == "Frontal Cortex")


set.seed(1234)
cere.dat <- cere.dat[order(cere.dat$Meth.Median),]
cere.mix <- normalmixEM(cere.dat$Meth.Median, mu=c(0,0.4, 0.6,1), k=4)
post.df <- as.data.frame(cbind(cere.mix$x,cere.mix$posterior))
colnames(post.df)[1] <- c("Methylation")
post.df <- post.df[order(post.df$Methylation),]
mid.comp <- post.df[((post.df$comp.2 > post.df$comp.1) &
                       (post.df$comp.2 > post.df$comp.3) &
                       (post.df$comp.2 > post.df$comp.4)),]
cere_methylation_cutoffs <- range(mid.comp$Methylation)


# Frontal
set.seed(1234)
frontal.dat <- frontal.dat[order(frontal.dat$Meth.Median),]
mix <- normalmixEM(frontal.dat$Meth.Median, mu=c(0,0.4, 0.6,1), k=4)
post.df <- as.data.frame(cbind(mix$x,mix$posterior))
colnames(post.df)[1] <- c("Methylation")
post.df <- post.df[order(post.df$Methylation),]
mid.comp <- post.df[((post.df$comp.2 > post.df$comp.1) &
                       (post.df$comp.2 > post.df$comp.3) &
                       (post.df$comp.2 > post.df$comp.4)),]
frontal_methylation_cutoffs <- range(mid.comp$Methylation)


################################
# Add the methylation clusters #
################################

dat$Methylation.Cluster <- "Medium Methylation"
dat$Methylation.Cluster[dat$Tissue == "Cerebellum" &
                          dat$Meth.Median < cere_methylation_cutoffs[1]] <- "Low Methylation"
dat$Methylation.Cluster[dat$Tissue == "Cerebellum" &
                          dat$Meth.Median > cere_methylation_cutoffs[2]] <- "High Methylation"

dat$Methylation.Cluster[dat$Tissue == "Frontal Cortex" &
                          dat$Meth.Median < frontal_methylation_cutoffs[1]] <- "Low Methylation"
dat$Methylation.Cluster[dat$Tissue == "Frontal Cortex" &
                          dat$Meth.Median > frontal_methylation_cutoffs[2]] <- "High Methylation"

write.table(dat, "Data/Outputs/Final Meth Summary.txt")



#######################
# Enrichment Analysis #
#######################

check_enrich <- function(dat) {
  df <- table(dat$Class, dat$Methylation.Cluster)
  resi <- list("Residuals" = round(chisq.test(df)$stdres,2),
               "P-Value" = chisq.test(df)$p.value)
  resi
}

meth_enrichments <- by(dat, INDICES = dat$Tissue, check_enrich)
meth_enrichments
