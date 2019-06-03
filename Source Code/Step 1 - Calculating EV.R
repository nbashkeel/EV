# This script transforms the raw datasets into an EV gene set classifications
setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(limma)
library(rowr)
library(plyr)
library(mixtools)
library(reshape2)
library(illuminaio)


####################################################################
# Section 1: Subset Brain Dataset into Cerebellum & Frontal Cortex #
####################################################################
brain.annotations <- read.csv("Data/Full Brain Annotations.csv")
colnames(brain.annotations) <- as.character(unlist(brain.annotations[2,]))
brain.annotations <- brain.annotations[c(1,2,10,11,12,13,14,15),-1]
brain.annotations <- as.data.frame(t(brain.annotations))
colnames(brain.annotations) <- c("Title", "Accession", "Tissue", "Sex", "Age", "PMI", "Tissuebank", "Batch")
samples <- rownames(brain.annotations)

brain.annotations <- apply(brain.annotations, 2, function(dat) {
  dat <- as.character(dat)
  dat <- gsub(".*: ","", dat)
})

brain.annotations <- as.data.frame(brain.annotations)
rownames(brain.annotations) <- samples

length(unique(paste0(brain.annotations$Tissuebank, brain.annotations$Batch)))
brain.annotations$Batch.Index <- as.factor(paste0(brain.annotations$Tissuebank, brain.annotations$Batch))

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


full.brain <- read.table("Data/brain_ds.txt")
full.brain <- full.brain[na.omit(match(rownames(genelist),rownames(full.brain))),]

brain.annotations <- brain.annotations[order(match(brain.annotations$Accession,colnames(full.brain))),]
batch.brain <- removeBatchEffect(full.brain, brain.annotations$Batch.Index)

cerebellum <- subset(brain.annotations, brain.annotations$Tissue == "cerebellum")
cerebellum <- batch.brain[,na.omit(match(rownames(cerebellum),colnames(batch.brain)))]
frontal <- subset(brain.annotations, brain.annotations$Tissue == "frontal cortex")
frontal <- batch.brain[,na.omit(match(rownames(frontal),colnames(batch.brain)))]

breast <- read.table("Data/normals_ExpressionMatrix.txt")
breast <- breast[na.omit(match(rownames(genelist),rownames(breast))),]


###################################################
# Section 2: Separate Bimodal & Unimodal Genes    #
#    End Result: Cerebellum Bimodal // Unimodal   #  
#                Frontal Bimodal // Unimodal      #  
#                Breast Bimodal // Unimodal       #  
###################################################


peaks.summary <- data.frame()
Bimodality <- function(dat){ 
  gene.name <- rownames(dat)
  dens <- na.omit(unlist(dat))
  dens <- dens[order(dens)]
  two.dist <- try(normalmixEM(dens,k=2,mu=c(min(dens),max(dens)),sigma=c(sd(dens),sd(dens)),arbvar=F,arbmean=F))
  if (two.dist == "Error in normalmixEM(dens, k = 2, mu = c(min(dens), max(dens)), sigma = c(sd(dens),  : \n  Too many tries!\n") {
    summary.peaks <- cbind(0,0,0,0,0,0,0,0,0)
    rownames(summary.peaks) <- gene.name
    colnames(summary.peaks) <- c("Peak 1 X", "Peak 2 X", "Peak 1 Y", "Peak 2 Y",
                                 "N1", "N2", "Peak Distance", "MAD", "Peak Ratio")
    peaks.summary <- rbind(peaks.summary,summary.peaks)}
  
  else {
    two.dist.post <- as.data.frame(cbind(two.dist$x, two.dist$posterior))
    two.dist.post$Diff <- abs(two.dist.post$comp.1-two.dist.post$comp.2)
    dens3 <- density(dens,n=length(dens))
    dens3 <- data.frame(X=dens3$x, Y=dens3$y)
    dens1 <- subset(dens3, dens3$X < two.dist.post[which.min(two.dist.post$Diff),1])
    dens2 <- subset(dens3, dens3$X > two.dist.post[which.min(two.dist.post$Diff),1])
    peaks <- c(dens1[which.max(dens1$Y),1], dens2[which.max(dens2$Y),1])
    
    peaks.Y <- c(dens3[which.min(abs(dens3$X-peaks[1])),2],
                 dens3[which.min(abs(dens3$X-peaks[2])),2]) 
    peak.distance <- abs(peaks[1]-peaks[2])
    peak.ratio <- min(peaks.Y)/max(peaks.Y)
    n1 <- length(subset(dens, dens < two.dist.post[which.min(two.dist.post$Diff),1]))
    n2 <- length(subset(dens, dens > two.dist.post[which.min(two.dist.post$Diff),1]))
    summary.peaks <- cbind(peaks[1], peaks[2], peaks.Y[1], peaks.Y[2],n1,n2,
                           peak.distance, mad(dens3$X),peak.ratio)
    rownames(summary.peaks) <- gene.name
    colnames(summary.peaks) <- c("Peak 1 X", "Peak 2 X", "Peak 1 Y", "Peak 2 Y", 
                                 "N1", "N2", "Peak Distance", "MAD", "Peak Ratio")
    peaks.summary <- rbind(peaks.summary,summary.peaks)
  }
}

Bimodal.DF <- function(dat) {
  test.dat <- apply(dat,1,Bimodality)
  test.dat <- do.call("rbind", test.dat)
  colnames(test.dat) <- c("Peak 1 X", "Peak 2 X", "Peak 1 Y", "Peak 2 Y",
                          "N1", "N2", "Peak Distance", "MAD", "Peak Ratio")
  test.dat$Bimodal <- NA
  test.dat$Bimodal[test.dat$`Peak Distance` > test.dat$MAD & 
                     test.dat$`Peak Ratio` > 0.1 & 
                     test.dat$N1 > (ncol(dat)/10) &
                     test.dat$N2 > (ncol(dat)/10)] <- "Y"
  test.dat
}

peaks.summary <- data.frame()
cere.bimodal <- Bimodal.DF(cerebellum)
cere.bimodal.y <- subset(cere.bimodal, cere.bimodal$Bimodal == "Y")
write.table(cere.bimodal.y, "Data/Outputs/Cerebellum Bimodal.txt")
cerebellum.bimodal <- cere.bimodal[!complete.cases(cere.bimodal),]
cerebellum <- as.data.frame(cerebellum[na.omit(match(rownames(cerebellum.bimodal),rownames(cerebellum))),])
write.table(cerebellum, "Data/Outputs/Cerebellum Unimodal.txt")

peaks.summary <- data.frame()
frontal.bimodal <- Bimodal.DF(frontal)
frontal.bimodal.y <- subset(frontal.bimodal, frontal.bimodal$Bimodal == "Y")
write.table(frontal.bimodal.y, "Data/Outputs/Frontal Bimodal.txt")
frontal.bimodal <- frontal.bimodal[!complete.cases(frontal.bimodal),]
frontal <- as.data.frame(frontal[na.omit(match(rownames(frontal.bimodal),rownames(frontal))),])
write.table(frontal, "Data/Outputs/Frontal Unimodal.txt")

peaks.summary <- data.frame()
breast.bimodal <- Bimodal.DF(breast)
breast.bimodal.y <- subset(breast.bimodal, breast.bimodal$Bimodal == "Y")
write.table(breast.bimodal.y, "Data/Outputs/Breast Bimodal.txt")
breast.bimodal <- breast.bimodal[!complete.cases(breast.bimodal),]
breast <- breast[na.omit(match(rownames(breast.bimodal),rownames(breast))),]
write.table(breast, "Data/Outputs/Breast Unimodal.txt")

################################################
# Section 3: Bootstrapping & EV Classification #
################################################
breast <- read.table("Data/Outputs/Breast Unimodal.txt")
cerebellum <- read.table("Data/Outputs/Cerebellum Unimodal.txt")
frontal <- read.table("Data/Outputs/Frontal Unimodal.txt")

Summarize.DF <- function(dat) {
  mad.dat <- apply(dat,1,function(dat) {
    boot <- numeric(1000)
    for (i in 1:1000) boot[i] <- mad(sample(as.vector(unlist(dat)),replace=T))
    mad.corrected <-  mad(dat) - (mean(boot)-mad(dat))
  })
  mad.dat <- data.frame(Obs.MAD=mad.dat)
  mad.dat$Probes <- rownames(mad.dat)
  
  dat$Probes <- rownames(dat)
  melted <- melt(dat, id.vars="Probes")
  dat <- ddply(melted, "Probes", summarise, Median=median(value,na.rm=T))
  dat <- merge(dat, mad.dat, by="Probes")
  dat <- dat[order(dat$Median),]
  dat
}


Quick_Summary <- function(dat){
  mad.dat <- apply(dat,1,function(dat) {
    gene.MAD <- mad(sample(as.vector(unlist(dat)),replace=T))
    mad.corrected <-  mad(dat) - (mean(gene.MAD)-mad(dat))
  })
  mad.dat <- data.frame(Obs.MAD=mad.dat)
  mad.dat$Probes <- rownames(mad.dat)
  
  dat$Probes <- rownames(dat)
  melted <- melt(dat, id.vars="Probes")
  dat <- ddply(melted, "Probes", summarise, Median=median(value,na.rm=T))
  dat <- merge(dat, mad.dat, by="Probes")
  dat <- dat[order(dat$Median),]
  dat
}

EV.dat <- function(dat){
  dat.exp <- loess(dat$Obs.MAD~dat$Median)
  dat.exp <- data.frame(dat.exp$x, dat.exp$fitted)
  colnames(dat.exp) <- c("Median", "Exp.MAD")
  dat <- merge(dat, dat.exp, by="Median")
  dat$EV <- dat$Obs.MAD-dat$Exp.MAD
  dat}

EVfunction <- function(dat2) {
  dat <- dat2$EV
  boot <- numeric(1000)
  for (i in 1:1000) boot[i] <- mad(sample(as.vector(unlist(dat)),replace=T))
  mad.corrected <-  mad(dat) - (mean(boot)-mad(dat))
  dat2.EV.range <- c(median(dat)-3*mad.corrected,
                     median(dat)+3*mad.corrected)
  dat2$Class <- "NV"
  dat2$Class[dat2$EV > dat2.EV.range[2]] <- "Hypervariable"
  dat2$Class[dat2$EV < dat2.EV.range[1]] <- "Hypovariable"
  dat2
}


Full.DF <- function(dat, Tissue) {
  dat <- Summarize.DF(dat)
  dat <- EV.dat(dat)
  dat <- EVfunction(dat)
  filestr <- paste0("Data/Outputs/", Tissue, " P1.txt")
  write.table(dat, filestr)
}

Full.DF(breast, "Breast")
Full.DF(cerebellum, "Cerebellum")
Full.DF(frontal, "Frontal")


#########################
# Section 4: Validation #
#########################

breast <- read.table("Data/Outputs/Breast Unimodal.txt")
cerebellum <- read.table("Data/Outputs/Cerebellum Unimodal.txt")
frontal <- read.table("Data/Outputs/Frontal Unimodal.txt")

breast.init <- data.frame("Probes" = rownames(breast))
cerebellum.init <- data.frame("Probes" = rownames(cerebellum))
frontal.init <- data.frame("Probes" = rownames(frontal))


Validation_function <- function(dat){
  mad.dat <- apply(dat,1,function(dat) {
    gene.MAD <- mad(sample(as.vector(unlist(dat)),replace=T))
    mad.corrected <-  mad(dat) - (mean(gene.MAD)-mad(dat))
  })
  mad.dat <- data.frame(Obs.MAD=mad.dat)
  mad.dat$Probes <- rownames(mad.dat)
  
  dat$Probes <- rownames(dat)
  melted <- melt(dat, id.vars="Probes")
  dat <- ddply(melted, "Probes", summarise, Median=median(value,na.rm=T))
  dat <- merge(dat, mad.dat, by="Probes")
  dat <- unique(dat)
  dat <- dat[order(dat$Median),]
  dat.exp <- loess(dat$Obs.MAD~dat$Median)
  dat.exp <- data.frame(dat.exp$x, dat.exp$fitted)
  colnames(dat.exp) <- c("Median", "Exp.MAD")
  dat <- merge(dat, dat.exp, by="Median")
  dat <- unique(dat)
  dat$EV <- dat$Obs.MAD-dat$Exp.MAD
  dat
  }

Validate <- function(dat, init){
  n <- ncol(dat)
  trainIndex <- sample(1:n, size = round(0.5*n), replace=FALSE)
  train <- dat[,trainIndex]
  train <- Validation_function(train)
  train <- unique(train)
  
  test <- dat[,-trainIndex]
  test <- Validation_function(test)
  test <- unique(test)
  
  train <- train[,c(2,5)]
  colnames(train)[2] <- c("Split1")
  test <- test[,c(2,5)]
  colnames(test)[2] <- c("Split2")
  
  df <- merge(train, test,by="Probes")
  return(df)
  }

for (i in 1:100){
  ptm <- proc.time()
  breast.init <- merge(breast.init, Validate(breast, breast.init), by="Probes")
  print(paste0("Breast Validation #",i,' Complete'))

  cerebellum.init <- merge(cerebellum.init, Validate(cerebellum, cerebellum.init), by="Probes")
  print(paste0("Cerebellum Validation #",i,' Complete'))
  
  frontal.init <- merge(frontal.init, Validate(frontal, frontal.init), by="Probes")
  print(paste0("Frontal Validation #",i,' Complete'))
  print(proc.time() - ptm)
}

write.table(breast.init, "Data/Validation/Breast.txt")
write.table(cerebellum.init, "Data/Validation/Cerebellum.txt")
write.table(frontal.init, "Data/Validation/Frontal.txt")


# Check how many hyper/hypo genes in each
class_table <- NULL
check_splits <- function(dat) {
  df <- data.frame("EV" = dat)
  df$Class <- "NV"
  df$Class[df$EV > (median(df$EV) + 3*mad(df$EV))] <- "Hypervariable"
  df$Class[df$EV < (median(df$EV) - 3*mad(df$EV))] <- "Hypovariable"
  class_table <- rbind(class_table, as.vector(table(df$Class)))
  class_table
}

breast.table <- apply(breast.init[,2:ncol(breast.init)], 2, check_splits)
cere.table <- apply(cerebellum.init[,2:ncol(cerebellum.init)], 2, check_splits)
frontal.table <- apply(frontal.init[,2:ncol(frontal.init)], 2, check_splits)

all.table <- as.data.frame(rbind(prop.table(rowMeans(breast.table))*100,
                                 prop.table(rowMeans(cere.table))*100,
                                 prop.table(rowMeans(frontal.table))*100))
colnames(all.table) <- c("Hypervariable", "Hypovariable", "NV")                   
rownames(all.table) <- c("Breast", "Cerebellum", "Frontal")

breast.class <- c(rep("Hypervariable", round(all.table[1,1]*nrow(breast.init)/100)),
                  rep("NV", round(all.table[1,3]*nrow(breast.init)/100)),
                  rep("Hypovariable", round(all.table[1,2]*nrow(breast.init)/100)))

cere.class <- c(rep("Hypervariable", round(all.table[2,1]*nrow(cerebellum.init)/100)),
                  rep("NV", round(all.table[2,3]*nrow(cerebellum.init)/100)),
                  rep("Hypovariable", round(all.table[2,2]*nrow(cerebellum.init)/100)))

frontal.class <- c(rep("Hypervariable", round(all.table[3,1]*nrow(frontal.init)/100)),
                rep("NV", round(all.table[3,3]*nrow(frontal.init)/100)-1),
                rep("Hypovariable", round(all.table[3,2]*nrow(frontal.init)/100)))

breast.init$Probes <- as.character(breast.init$Probes)
cerebellum.init$Probes <- as.character(cerebellum.init$Probes)
frontal.init$Probes <- as.character(frontal.init$Probes)


breast_df <- apply(breast.init[,2:ncol(breast.init)], 2, function(x) {
  test <- data.frame("Data" = x, "Probes" = breast.init$Probes)
  test <- test[order(-test$Data),]
  test$Class <- breast.class
  test <- test[order(test$Probes),]
  test$Class
})
genes <- breast.init[,'Probes',drop=F]
genes <- genes[order(genes$Probes), ,drop=F]
breast_df <- cbind(genes, breast_df)

cerebellum_df <- apply(cerebellum.init[,2:ncol(breast.init)], 2, function(x) {
  test <- data.frame("Data" = x, "Probes" = as.character(cerebellum.init$Probes))
  test <- test[order(-test$Data),]
  test$Class <- cere.class
  test <- test[order(test$Probes),]
  test$Class
})
genes <- cerebellum.init[,'Probes',drop=F]
genes <- genes[order(genes$Probes), ,drop=F]
cerebellum_df <- cbind(genes, cerebellum_df)

gene_class3 <- NULL
frontal_df <- apply(frontal.init[,2:ncol(breast.init)], 2, function(x) {
  test <- data.frame("Data" = x, "Probes" = as.character(frontal.init$Probes))
  test <- test[order(-test$Data),]
  test$Class <- frontal.class
  test <- test[order(test$Probes),]
  test$Class
})
genes <- frontal.init[,'Probes',drop=F]
genes <- genes[order(genes$Probes), ,drop=F]
frontal_df <- cbind(genes, frontal_df)

# Check if splits match & calculate gene mean
splits_mean <- function(dat, df){
  df <- df[,c(2,6)]
  colnames(df)[2] <- c("Original.EV")
  df <- df[order(df$Probes),]
  runs <- seq(2,ncol(breast_df),2)
  
  for (i in 1:length(runs)){
    df2 <- cbind(df[,1:2], dat[,runs[i]], dat[,runs[i]+1])
    df <- cbind(df, ((df2[,2] == df2[,3]) & (df2[,2] == df2[,4]))*1)
  }
  df
}


rowwise_binom <- function(x) {
  x = unname(unlist(x))
  binom.test(c(length(x[x == 0]), length(x[x == 1])), n=ncol(matching_df), 
             p=0.5, alternative="less")$p.value
}

breast.og <- read.table("Data/Outputs/Breast P1.txt")
breast.means <- splits_mean(breast_df, breast.og)
rownames(breast.means) <- breast.means$Probes
colnames(breast.means) <- c("Probes", "Original.EV", seq(1,ncol(breast.means)-2,1))
test.pvals <- apply(breast.means[,3:ncol(breast.means)], 1, rowwise_binom)
test.pvals <- data.frame("pvals" = test.pvals, "Probes" = names(test.pvals))
breast.means <- merge(test.pvals, breast.means, by="Probes")
breast.means$Average <- rowMeans(breast.means[,4:ncol(breast.means)])
breast.means <- breast.means[,c(1,2,3,ncol(breast.means))]
breast.means$Tissue <- "Breast"

cere.og <- read.table("Data/Outputs/Cerebellum P1.txt")
cere.means <- splits_mean(cerebellum_df, cere.og)
rownames(cere.means) <- cere.means$Probes
colnames(cere.means) <- c("Probes", "Original.EV", seq(1,ncol(cere.means)-2,1))
test.pvals <- apply(cere.means[,3:ncol(cere.means)], 1, rowwise_binom)
test.pvals <- data.frame("pvals" = test.pvals, "Probes" = names(test.pvals))
cere.means <- merge(test.pvals, cere.means, by="Probes")
cere.means$Average <- rowMeans(cere.means[,4:ncol(cere.means)])
cere.means <- cere.means[,c(1,2,3,ncol(cere.means))]
cere.means$Tissue <- "Cerebellum"

frontal.og <- read.table("Data/Outputs/Frontal P1.txt")
frontal.means <-splits_mean(frontal_df, frontal.og)
rownames(frontal.means) <- frontal.means$Probes
colnames(frontal.means) <- c("Probes", "Original.EV", seq(1,ncol(frontal.means)-2,1))
test.pvals <- apply(frontal.means[,3:ncol(frontal.means)], 1, rowwise_binom)
test.pvals <- data.frame("pvals" = test.pvals, "Probes" = names(test.pvals))
frontal.means <- merge(test.pvals, frontal.means, by="Probes")
frontal.means$Average <- rowMeans(frontal.means[,4:ncol(frontal.means)])
frontal.means <- frontal.means[,c(1,2,3,ncol(frontal.means))]
frontal.means$Tissue <- "Frontal Cortex"

df <- rbind(breast.means, cere.means, frontal.means)
df$Sig <- "NS"
df$Sig[df$pvals < 0.05] <- "Sig"

colnames(df) <- c("Probes", "pvals", "Class", "Means", "Tissue", "Sig")
genelist <- readBGX("Data/HumanHT-12_V3_0_R3_11283641_A.bgx")
genelist <- data.frame(genelist$probes)
genelist <- data.frame(Gene.Symbol=genelist$ILMN_Gene,
                       Probes=genelist$Probe_Id)
df <- merge(genelist, df, by="Probes")
write.table(df, "Data/Outputs/Final All 3.txt")