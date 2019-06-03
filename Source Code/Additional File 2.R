setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))
library(ggplot2)
library(grid)
library(gridExtra)
library(ggplot2)
library(illuminaio)
library(rowr)
library(plyr)
library(reshape2)
library(lsa)
library(gtools)
library(viridis)


df <- read.table("Data/Outputs/Final All 3.txt")
breast <- read.table("Data/Outputs/Breast P1.txt")
breast$Tissue <- "Breast"
cerebellum <- read.table("Data/Outputs/Cerebellum P1.txt")
cerebellum$Tissue <- "Cerebellum"
frontal <- read.table("Data/Outputs/Frontal P1.txt")
frontal$Tissue <- "Frontal Cortex"
EVs <- rbind(breast, cerebellum, frontal)
EVs <- EVs[,c(2,5,7)]


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
genelist$Coordinates <- gsub("-",",",genelist$Coordinates)
genelist$Coordinates <- gsub(":",",",genelist$Coordinates)

for (i in 1:nrow(genelist)){
  texta <- paste0("c(",genelist[i,6],")")
  texta <- eval(parse(text = texta))
  texta <- round(median(texta),0)
  genelist[i,6] <- texta
} # Max range is 999 bp, basically negligible

genelist <- unique(genelist[,c(1,5,6)]) # 91 probes removed (41993 probes)
genelist$Chromosome <- as.character(genelist$Chromosome)
genelist <- subset(genelist, nchar(genelist$Chromosome) < 4) # 107 probes removed (41886 probes)
genelist <- subset(genelist, genelist$Chromosome != "XY" & 
                     genelist$Chromosome != "YX") # 9 Probes removed
genelist$Probes <- rownames(genelist)
genelist$Coordinates <- as.numeric(genelist$Coordinates)
genelist <- genelist[with(genelist, order(Chromosome,Coordinates)),]
genelist$Gene.Symbol <- as.character(genelist$Gene.Symbol)
x <- rle(genelist$Gene.Symbol)
y <- cumsum(c(1, x$lengths[-length(x$lengths)]))
genelist <- genelist[y,]
chrs <- data.frame(Max=tapply(genelist$Coordinates, genelist$Chromosome, max))
chrs$Chromosome <- as.factor(rownames(chrs))

genelist <- merge(genelist, chrs, by="Chromosome")

genelist$Max <- as.vector(genelist$Max)
genelist$Bin <- ceiling(genelist$Coordinates/(genelist$Max/100))
genelist$Bin[genelist$Bin > 100] <- 100
genelist$Chr.Mod <- as.character(genelist$Chromosome)
genelist$Chr.Mod <- revalue(genelist$Chr.Mod, replace=c("X" = 23, "Y" = 24))
genelist$Chr.Mod <- as.numeric(genelist$Chr.Mod)
genelist$Bin.Mod <- genelist$Bin+((genelist$Chr.Mod-1)*100)


df <- merge(df, EVs, by=c("Probes","Tissue"))
df <- subset(df, df$Sig == "Sig")

DF.summary <- function(dat) {
  dat <- merge(dat, genelist[,c(2,8)], by="Gene.Symbol")
  DF.Summary <- ddply(dat, ~Bin.Mod, summarise, EV.Bin=mean(EV))
  DF.Summary
}

breast <- DF.summary(subset(df, df$Tissue == "Breast"))
breast$Tissue <- "Breast"
cerebellum <- DF.summary(subset(df, df$Tissue == "Cerebellum"))
cerebellum$Tissue <- "Cerebellum"
frontal <- DF.summary(subset(df, df$Tissue == "Frontal Cortex"))
frontal$Tissue <- "Frontal"

ggDF <- rbind(breast, cerebellum, frontal)
ggDF$Sign <- as.factor(sign(ggDF$EV.Bin))
ranges <- data.frame("chrs" = seq(1,2400, by=100))
ranges2 <- seq(51,2450, by=100)
xlabels <- c(1:22, "X", "Y")


p3 <- ggplot(ggDF,aes(x=Bin.Mod,xend=Bin.Mod,y=0,yend=EV.Bin, col=Sign))+
  geom_segment(size=0.25) + 
  facet_wrap(~Tissue, scales="free", nrow=3, strip.position = "right") +
  xlab("Chromosome") + ylab("EV")+
  theme(axis.text.y=element_text(size=8),
        axis.title.x = element_text(vjust=5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_continuous(limits=c(0,2400), expand = c(0, 0),
                     breaks = ranges2, labels=xlabels) +
  geom_hline(aes(yintercept=0), alpha=0.25, lwd=0.25, lty=5) +
  geom_vline(data=ranges, aes(xintercept=chrs), alpha=0.5, lty=3) +
  scale_color_grey(start=0.5, end=0) +
  guides(color=FALSE)

p3





pos.breast <- subset(ggDF, ggDF$Tissue == "Breast")
#pos.breast$EV.Bin[pos.breast$EV.Bin < 0] <- 0
pos.breast <- pos.breast[,-4]

chr_df <- data.frame("Chr" = rep(c(1:22, "X", "Y"), each = 100))
chr_df$Bin.Mod <- seq(1,2400,1)

pos.breast <- merge(pos.breast, chr_df, by="Bin.Mod", all.y = TRUE)
pos.breast$EV.Bin[is.na(pos.breast$EV.Bin)] <- 0
pos.breast$Tissue[is.na(pos.breast$Tissue)] <- 'Breast'
pos.breast$Bin.Mod <- rep(seq(1,100,1),24)

recast.pos.breast <- dcast(data = pos.breast, formula = pos.breast$Bin.Mod~pos.breast$Chr, value.var = "EV.Bin")
recast.pos.breast <- as.matrix(recast.pos.breast[,-1])

test <- cosine(recast.pos.breast)
v.exclude <- c(0,1)
max.exclude <- max(test[!test %in% v.exclude])
min.exclude <- min(test[!test %in% v.exclude])
c(min.exclude, max.exclude)


clean_chr_df <- function(tissue){
  tissue_df <- subset(ggDF, ggDF$Tissue == tissue)
  #tissue_df$EV.Bin[tissue_df$EV.Bin < 0] <- 0
  tissue_df <- tissue_df[,-4]
  
  chr_df <- data.frame("Chr" = rep(c(1:22, "X", "Y"), each = 100))
  chr_df$Bin.Mod <- seq(1,2400,1)
  
  tissue_df <- merge(tissue_df, chr_df, by="Bin.Mod", all.y = TRUE)
  tissue_df$EV.Bin[is.na(tissue_df$EV.Bin)] <- 0
  tissue_df$Tissue[is.na(tissue_df$Tissue)] <- 'Breast'
  tissue_df$Bin.Mod <- rep(seq(1,100,1),24)
  
  recast.tissue_df <- dcast(data = tissue_df, formula = tissue_df$Bin.Mod~tissue_df$Chr, value.var = "EV.Bin")
  recast.tissue_df <- as.matrix(recast.tissue_df[,-1])
  recast.tissue_df <- recast.tissue_df[,mixedorder(colnames(recast.tissue_df))]
  recast.tissue_df
  }




clean_breast <- clean_chr_df("Breast")
clean_cerebellum <- clean_chr_df("Cerebellum")
clean_frontal <- clean_chr_df("Frontal")

library(gplots)

similarity_matrix <- function(dat, tissue){
  test <- cosine(dat)
  test <- test[mixedorder(rownames(test)),mixedorder(colnames(test))]
  v.exclude <- c(0,1)
  max.exclude <- max(test[!test %in% v.exclude])
  min.exclude <- min(test[!test %in% v.exclude])
  print(c(min.exclude, max.exclude))
  my_palette <- colorRampPalette(c("red", "black", "green"))(n = 300)
  test[test == 1] <- 0
  title_string <- paste0("Chromosome Correlation\nin ", toString(tissue), " Tissue\n",
                         "Min: ", format(min.exclude,digits=3), "  Max: ", format(max.exclude,digits=3))
  
  heatmap.2(test,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
            col = my_palette, density.info ='none', srtCol=0, main = title_string)
}

similarity_matrix(clean_breast, "Breast")
similarity_matrix(clean_cerebellum, "Cerebellum")
similarity_matrix(clean_frontal, 'Frontal')




tissue_comparisons <- NULL
for (i in 1:24){
  all_tissues <- cbind(clean_breast[,i], clean_cerebellum[,i], clean_frontal[,i])
  colnames(all_tissues) <- c("Breast", "Cerebellum", "Frontal")
  all_cosine <- cosine(all_tissues)
  tissue_comparisons <- rbind(tissue_comparisons,all_cosine[lower.tri(all_cosine)])
}

colnames(tissue_comparisons) <- c("Breast-Cerebellum", "Breast-Frontal", "Cerebellum-Frontal")
rownames(tissue_comparisons) <- colnames(clean_breast)

heatmap.2(tissue_comparisons,symkey=FALSE, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          col = viridis(100), density.info ='none', srtCol=0, main="Comparing Chromosomes\nAcross Tissue Types",
          cexCol = 1.25, adjCol = 0.5, symbreaks = FALSE)



pdf("Figures/Additional File 2.pdf")
textplot("Additional File 2. EV correlation between different tissue types
", halign="center", valign='top')
similarity_matrix(clean_breast, "Breast")
similarity_matrix(clean_cerebellum, "Cerebellum")
similarity_matrix(clean_frontal, 'Frontal')

heatmap.2(tissue_comparisons,symkey=FALSE, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          col = viridis(100), density.info ='none', srtCol=0, main="Comparing Chromosomes\nAcross Tissue Types",
          cexCol = 1.25, adjCol = 0.5, symbreaks = FALSE)
dev.off()
