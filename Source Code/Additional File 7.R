setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(illuminaio)
library(pvclust)
library(mclust)
library(gplots)

################
# Age Clusters #
################

age <- read.csv("Data/Full Brain Annotations.csv")
colnames(age) <- as.vector(unlist(age[2,]))
age <- age[,-1]
age <- age[c(10,12),]
age <- as.data.frame(t(age))
samples <- rownames(age)

age <- apply(age, 2, function(dat) {
  dat <- as.character(dat)
  dat <- gsub(".*: ", "", dat)
})

age <- as.data.frame(age)
rownames(age) <- samples
colnames(age) <- c("Tissue", "Age")
age$Age <- as.numeric(age$Age)

set.seed(123)
k.max <- 15
scaled_age <- as.matrix(scale(age$Age))
age_all <- sapply(1:k.max, function(k){kmeans(scaled_age, k, nstart=50,iter.max = 15 )$tot.withinss})



#################################################################
# Hierarchical Clustering of Gene Expression Brain Data Subsets #
#################################################################

illumina <- readBGX("Data/HumanHT-12_V3_0_R3_11283641_A.bgx")
illumina_probes <- data.frame(illumina$probes)
genelist <- data.frame(Gene.Symbol=illumina_probes$ILMN_Gene,
                       Probes=illumina_probes$Probe_Id,
                       Refseq=illumina_probes$RefSeq_ID,
                       Unigene=illumina_probes$Unigene_ID)
genelist <- genelist[-grep("NR",genelist$Refseq),] # 654 genes
genelist <- genelist[-grep("XR",genelist$Refseq),] # 501 genes
rownames(genelist) <- genelist$Probe
genelist$Probes <- rownames(genelist)

test <- read.csv("Data/Full Brain Annotations.csv")
colnames(test) <- as.vector(unlist(test[2,]))
test <- test[,-1]
test <- test[c(10,11,12,13),]
test <- as.data.frame(t(test))
colnames(test) <- c("Tissue", "Sex", "Age", "PMI")
samples <- rownames(test)

test <- apply(test, 2, function(dat) {
  dat <- as.character(dat)
  dat <- gsub(".*: ", "", dat)
})

test <- as.data.frame(test)
rownames(test) <- samples
test$Age <- as.numeric(test$Age)
test$Age.Cat <- 1
test[test$Age < 22,5] <- 1
test[test$Age >= 22 & test$Age < 74,5] <- 2
test[test$Age >= 74 ,5] <- 3
test$Age.Cat <- as.factor(test$Age.Cat)
test$Age <- as.numeric(test$Age)
test <- test[order(rownames(test)),]
test$PMI <- as.numeric(test$PMI)

uniquerows <- test[,c(1,2,5)]
uniquerows <- unique(uniquerows)

quicknames <- uniquerows
quicknames$Tissue <- gsub("cerebellum", "C", quicknames$Tissue)
quicknames$Tissue <- gsub("frontal cortex", "FC", quicknames$Tissue)
quicknames$Sex <- gsub("male","M", quicknames$Sex)
quicknames$Sex <- gsub("feM","F", quicknames$Sex)
for (i in 1:nrow(quicknames)) {
  rownames(quicknames)[i] <- paste0(quicknames[i,1],quicknames[i,2],quicknames[i,3])
  rownames(uniquerows)[i] <- paste0(quicknames[i,1],quicknames[i,2],quicknames[i,3])
}

full.brain <- read.table("Data/brain_ds.txt")
full.brain <- full.brain[na.omit(match(rownames(genelist),rownames(full.brain))),]

EV.df <- as.data.frame(matrix(NA,nrow=nrow(full.brain),ncol=0))
EV.df$Probes <- rownames(full.brain)

for (i in 1:nrow(uniquerows)) {
  tissue <- uniquerows[i,1]
  sex <- uniquerows[i,2]
  age <- uniquerows[i,3]
  
  subdat <- subset(test, test$Tissue == tissue & 
                     test$Sex == sex & test$Age.Cat == age)
  brain.dat <- full.brain[,na.omit(match(rownames(subdat),colnames(full.brain)))]
  brain.dat <- data.frame(apply(brain.dat,1,median))
  brain.dat$Probes <- rownames(brain.dat)
  EV.df <- merge(EV.df, brain.dat, by="Probes")
  col <- i+1
  colnames(EV.df)[col] <- paste0(quicknames[i,1],quicknames[i,2],quicknames[i,3])
}
rownames(EV.df) <- EV.df$Probes
EV.df <- EV.df[,-1]


fit <- pvclust(EV.df, method.hclust="ward",
               method.dist="euclidean")



pdf("Figures/Additional File 7.pdf")
textplot("Additional File 7. Preprocessing of Brain Samples

1) K-Means clustering of brain sample age to determine 
  the optimal number of age clusters using the within-clusters
  sum of squares elbow method
  
  
2) Hierarchical clustering via multiscale bootstrap of 
  12 brain sample subset permutations using pvclust package
    12 Groups = Age Category(3) x Sex(2) x Tissue Type(2)
", halign="center", valign='top')

plot(1:k.max, age_all,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main = "Kmeans Clustering of Age of Brain Samples")
abline(v=3, col='red', lty=2)
plot(fit, print.num = FALSE, float=0.03)
pvrect(fit, alpha=0.99)
dev.off()
