setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(illuminaio)
library(grid)
library(gridExtra)
library(gplots)



df <- read.table("Data/Outputs/Final All 3.txt")
df$Class[df$Sig == "NS"] <- "NV"
breast <- subset(df, df$Tissue == "Breast")
cerebellum <- subset(df, df$Tissue == "Cerebellum")
frontal <- subset(df, df$Tissue == "Frontal Cortex")


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

genelist <- unique(genelist[,c(1,2,5,6)]) # 91 probes removed (41993 probes)
genelist$Chromosome <- as.character(genelist$Chromosome)
genelist <- subset(genelist, nchar(genelist$Chromosome) < 4) # 107 probes removed (41886 probes)
genelist <- subset(genelist, genelist$Chromosome != "XY" & 
                     genelist$Chromosome != "YX") # 9 Probes removed



Essentiality <- function(dat, table_title) {
  essentialgenes <- read.csv("Data/Essential Genes.csv")
  essentialgenes <- subset(essentialgenes, essentialgenes$Essential == "Y" & essentialgenes$Non.essential == "N")
  essentialgenes <- essentialgenes[,c(1,6)]
  colnames(essentialgenes) <- c("Gene.Symbol", "Essential")
  
  nr <- nrow(dat)
  dat$Essential <- "Non-Essential"
  dat[na.omit(match(essentialgenes$Gene.Symbol, dat$Gene.Symbol)),"Essential"] <- "Essential"
  dat$Essential <- as.factor(dat$Essential)
  
  counts <- table(dat$Essential, dat$Class)
  stdres <- format(chisq.test(counts)$stdres, digits=4)
  pval <- chisq.test(counts)$p.value
  chi_summary <- chisq.test(counts)
  chi_summary <- data.frame(chi_summary[[1]],chi_summary[[2]], chi_summary[[3]])
  colnames(chi_summary) <- c("X-Squared", "Degrees of Freedom", "P-Value")
  rownames(chi_summary) <- c("")
  
  counts <- rbind(counts,colSums(counts))
  counts <- cbind(counts, rowSums(counts))
  
  
  grob(textplot(counts))
  title(paste("\n", table_title, "Tissue Essential Gene Enrichment Analysis\nContingency Table"))
  grob(textplot(chi_summary))
  title("Chi-Squared Test Summary")
  grob(textplot(stdres))
  title("Standardized Residuals")
}

###############
# Methylation #
###############

meth_dat<- read.table("Data/Outputs/Final Meth Summary.txt")
meth_cere <- subset(meth_dat, meth_dat$Tissue == "Cerebellum")
meth_fc <- subset(meth_dat, meth_dat$Tissue == "Frontal Cortex")

Methylation_chi <- function(dat, table_title) {
  counts <- table(dat$Class, dat$Methylation.Cluster)
  stdres <- format(chisq.test(counts)$stdres, digits=4)
  pval <- chisq.test(counts)$p.value
  chi_summary <- chisq.test(counts)
  chi_summary <- data.frame(chi_summary[[1]],chi_summary[[2]], chi_summary[[3]])
  colnames(chi_summary) <- c("X-Squared", "Degrees of Freedom", "P-Value")
  rownames(chi_summary) <- c("")
  
  counts <- rbind(counts,colSums(counts))
  counts <- cbind(counts, rowSums(counts))
  
  grob(textplot(counts))
  title(paste("\n", table_title, "Tissue Methylation Enrichment Analysis\nContingency Table"))
  grob(textplot(chi_summary))
  title("Chi-Squared Test Summary")
  grob(textplot(stdres))
  title("Standardized Residuals")
}

Methylation_chi(meth_cere, "Cerebellum")
Methylation_chi(meth_fc, "Frontal Cortex")



pdf("Figures/Additional File 4.pdf")
par(mfrow=c(2,1))
textplot("Additional File 4. Chi-Squared Enrichment Analysis\nMethodology
", halign="center", valign='top')
textplot("
  1) Construct a contigency table of genes for each tissue type
  2) Conduct a chi-square test using the chisq.test() function in R
  3) Extract the standardized residuals(stdres) whereby:\n
  > standardized residuals, (observed - expected) / sqrt(V), where V is the
  > residual cell variance (Agresti, 2007, section 2.4.5 for the case where 
  > x is a matrix, n * p * (1 - p) otherwise).", cex=0.75, valign="top")


par(mfrow=c(3,1))
Essentiality(breast, "Breast")
Essentiality(cerebellum, "Cerebellum")
Essentiality(frontal, "Frontal Cortex")

Methylation_chi(meth_cere, "Cerebellum")
Methylation_chi(meth_fc, "Frontal Cortex")
dev.off()
