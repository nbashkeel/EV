setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))
library(plyr)
library(plotrix)
library(gplots)
library(illuminaio)

df <- read.table("Data/Outputs/Final All 3.txt")
df <- subset(df, df$Sig == "Sig" & df$Class == "Hypervariable")

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
test$PMI <- as.numeric(test$PMI)
test$Sex <- as.factor(test$Sex)

age.cere <- subset(test, test$Tissue == "cerebellum")
cerebellum <- subset(df, df$Tissue == "Cerebellum")
cere.hyper <- read.table("Data/Outputs/Cerebellum Unimodal.txt")
cere.hyper <- cere.hyper[na.omit(match(cerebellum$Probes,rownames(cere.hyper))),]

age.frontal <- subset(test, test$Tissue == "frontal cortex")
frontal <- subset(df, df$Tissue == "Frontal Cortex")
frontal.hyper <- read.table("Data/Outputs/Frontal Unimodal.txt")
frontal.hyper <- frontal.hyper[na.omit(match(frontal$Probes, rownames(frontal.hyper))),]

genelist <- readBGX("Data/HumanHT-12_V3_0_R3_11283641_A.bgx")
genelist <- data.frame(genelist$probes)
genelist <- data.frame(Gene.Symbol=genelist$ILMN_Gene,
                       Probes=genelist$Probe_Id)

dat <- NULL
colIndex <- NULL
draw.heatmap <- function(DF, direction){
  if (DF == "cere") {dat <- cere.hyper}
  else if (DF == "frontal") {dat <- frontal.hyper}
  else {print ("Incorrect dataframe")}
  
  dat <- t(dat)
  dat <- merge(test, dat, by="row.names")
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  options(contrasts = c("contr.sum","contr.poly"))
  test.DF <- NULL
  for (i in 4:ncol(dat)) {
    gene <- colnames(dat)[i]
    test.lm <- summary(lm(dat[,i]~dat$PMI+dat$Sex+dat$Age))
    test.lm <- as.data.frame(cbind(test.lm$coefficients['dat$Age',4],
                                   sign(test.lm$coefficients['dat$Age',1]),
                                   gene))
    test.DF <- rbind(test.DF, test.lm)}
  colnames(test.DF) <- c("Age.Pr", "Age.Sign", "Probe")
  test.DF$Age.Pr <- as.numeric(as.character(test.DF$Age.Pr))
  test.DF$Age.Pr <- p.adjust(test.DF[order(test.DF$Age.Pr),"Age.Pr"], method="BH", n=nrow(test.DF))
  
  
  if (direction == "up") {
    test.DF <- subset(test.DF, test.DF$Age.Pr < 0.01 & 
                        test.DF$Age.Sign == 1)
  }
  else if (direction == "down") {
    test.DF <- subset(test.DF, test.DF$Age.Pr < 0.01 & 
                        test.DF$Age.Sign == -1)
  }
  else {print ("Incorrect Direction")}
  
  dat <- dat[,c(na.omit(match(test.DF$Probe, colnames(dat))))]
  set.seed(1234)
  dat <- scale(dat)
  dat[dat > 3] <- 3
  dat[dat < -3] <- -3
  dat <- t(dat)
  
  if (DF == "cere") {age.dist <- dist(matrix(age.cere$Age))}
  else if (DF == "frontal") {age.dist <- dist(matrix(age.frontal$Age))}
  else {print ("Incorrect dataframe")}
  
  age.clust <- hclust(age.dist)
  age.dend <- as.dendrogram(age.clust)
  age.clustcut <- cutree(age.clust, h=45)
  age.cols <- rainbow(3)
  age.cluster.cols <- age.cols[age.clustcut]
  
  colIndex <- table(age.clustcut)
  col.cols <- lapply(rainbow(length(unique(age.clustcut))),color.id)
  col.cols <- do.call(rbind,col.cols)
  names(colIndex) <- col.cols[,1]
  colIndex <- colIndex[c("blue", "green", "red")]
  
  gene.dist <- dist(as.matrix(dat))
  gene.clust <- hclust(gene.dist)
  gene.dend <- as.dendrogram(gene.clust)
  plot(gene.dend)
  repeat{
    n <- as.integer(readline(prompt="Cut gene dendrogram at: "))
    if (n == 0) {break}
    gene.clustcut <- cutree(gene.clust, h=n)
    gene.cols <- rainbow(15)[1:length(unique(gene.clustcut))]
    gene.cluster.cols <- gene.cols[gene.clustcut]
    
    h <- heatmap.2(dat,col=redgreen(100),
                   Colv=as.dendrogram(age.clust), Rowv=as.dendrogram(gene.clust),
                   ColSideColors = age.cluster.cols, RowSideColors = gene.cluster.cols,
                   labRow="", labCol="",
                   symm=F,symkey=F,symbreaks=F,scale="none", density.info="none", trace="none")
  }
  
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
  }
  
  if (DF == "cere") {DF2 <- "Cerebellum"}
  else if (DF == "frontal") {DF2 <- "Frontal Cortex"}
  
  filename <- paste0("Figures/",DF2, " ", simpleCap(direction),"regulated Age Genes.tiff")
  tiff(filename=filename, width=2250, height=(2250/2),
      units="px", pointsize=12)
  heatmap.2(dat,col=redgreen(100),
            Colv=as.dendrogram(age.clust), Rowv=as.dendrogram(gene.clust),
            ColSideColors = age.cluster.cols, RowSideColors = gene.cluster.cols,
            labRow="", labCol="",
            symm=F,symkey=F,symbreaks=F,scale="none", density.info="none", trace="none",
            add.expr = abline(v = c(155,313,456), lwd = 0.5, col="white"),
            key=F)
  dev.off()
  
  test.row <- gene.cluster.cols[h$rowInd]
  test.row <- sapply(test.row,color.id)
  test.row <- do.call(rbind,test.row)
  test.row <- data.frame(Color=test.row[,1], Probes=rownames(dat))
  test.row <- merge(test.row, genelist,by="Probes")
  
  for (i in 1:length(unique(test.row$Color))){
    color <- as.vector(unique(test.row$Color)[i])
    test2 <- subset(test.row, test.row$Color==color)
    file.cluster <- paste0("Data/Outputs/GO/Age/", DF2, " ", simpleCap(direction),"regulated Age Gene ",simpleCap(color), " Cluster.txt")
    write.table(test2$Gene.Symbol, file=file.cluster, quote=F, row.names = F, col.names = F)
  }
  
  col.summary <- as.data.frame(table(test.row$Color))
  colnames(col.summary) <- c("Color", "Frequency")
  col.summary <- col.summary[order(-col.summary$Frequency),]
  col.sum.name <- paste0("Data/Outputs/GO/Age/", DF2, " ", simpleCap(direction),"regulated Age Gene Color Summary.csv")
  write.csv(col.summary, file=col.sum.name)
}


draw.heatmap(DF="cere", direction="up") #Type in (33) when prompted, followed by (0) at second prompt
draw.heatmap(DF="cere", direction="down") #Type in (32) when prompted, followed by (0) at second prompt

draw.heatmap(DF="frontal", direction="up") #Type in (34) when prompted, followed by (0) at second prompt
draw.heatmap(DF="frontal", direction="down") #Type in (34) when prompted, followed by (0) at second prompt

# Merge the two plots in Photoshop