setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))
library(illuminaio)

df <- read.table("Data/Outputs/Final All 3.txt")
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



Essentiality <- function(dat) {
  essentialgenes <- read.csv("Data/Essential Genes.csv")
  essentialgenes <- subset(essentialgenes, essentialgenes$Essential == "Y" & essentialgenes$Non.essential == "N")
  essentialgenes <- essentialgenes[,c(1,6)]
  colnames(essentialgenes) <- c("Gene.Symbol", "Essential")
  
  dat <- subset(dat, dat$Sig == "Sig")
  nr <- nrow(dat)
  dat$Essential <- "N"
  dat[na.omit(match(essentialgenes$Gene.Symbol, dat$Gene.Symbol)),"Essential"] <- "Y"
  dat$Essential <- as.factor(dat$Essential)
  
  counts <- table(dat$Essential, dat$Class)
  stdres <- chisq.test(counts)$stdres
  pval <- chisq.test(counts)$p.value
  counts <- rbind(counts,colSums(counts))
  
  
  summary.table <- as.data.frame(matrix(c(counts[3,1], counts[2,1], round(stdres[2,1],2),
                                          counts[3,2], counts[2,2], round(stdres[2,2],2),
                                          counts[3,3], counts[2,3], round(stdres[2,3],2)),
                                        ncol=3, byrow=T))
  colnames(summary.table) <- c("N" , "Essential Gene Counts", "Standard Residuals")
  rownames(summary.table) <- c("Hyper", "Hypo", "NV")
  summary.table$`P Value` <- pval
  summary.table
}


breast.essential <- Essentiality(breast)
cere.essential <- Essentiality(cerebellum)
frontal.essential <- Essentiality(frontal)
