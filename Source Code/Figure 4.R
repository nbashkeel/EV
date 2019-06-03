setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))


library(ggplot2)
library(grid)
library(gridExtra)
library(VennDiagram)
library(illuminaio)
library(rowr)
library(plyr)
library(viridis)

df <- read.table("Data/Outputs/Final All 3.txt")
df$Class <- revalue(df$Class, c("NV" = "Non-Variable", "Hypervariable" = "Hyper-Variable",
                                "Hypovariable" = "Hypo-Variable"))


breast <- subset(df, df$Tissue == "Breast")
cerebellum <- subset(df, df$Tissue == "Cerebellum")
frontal <- subset(df, df$Tissue == "Frontal Cortex")

breast.hyper <- subset(breast, breast$Class=="Hyper-Variable" & breast$Sig == "Sig")
cerebellum.hyper <- subset(cerebellum, cerebellum$Class=="Hyper-Variable" & cerebellum$Sig == "Sig")
frontal.hyper <- subset(frontal, frontal$Class=="Hyper-Variable" & frontal$Sig == "Sig")

breast.hypo <- subset(breast, breast$Class=="Hypo-Variable" & breast$Sig == "Sig")
cerebellum.hypo <- subset(cerebellum, cerebellum$Class=="Hypo-Variable" & cerebellum$Sig == "Sig")
frontal.hypo <- subset(frontal, frontal$Class=="Hypo-Variable" & frontal$Sig == "Sig")

breast.NV <- subset(breast, breast$Class=="Non-Variable" | breast$Sig == "NS")
cerebellum.NV <- subset(cerebellum, cerebellum$Class=="Non-Variable" | cerebellum$Sig == "NS")
frontal.NV <- subset(frontal, frontal$Class=="Non-Variable" | frontal$Sig == "NS")


venn1 <- venn.diagram(x=list(
  "Breast" = breast.hyper$Probes,
  "Cerebellum" = cerebellum.hyper$Probes,
  "Frontal Cortex" = frontal.hyper$Probes), file=NULL,
  scaled = FALSE,fill = viridis(3), main="Hypervariable Gene Sets", cat.dist=0.09,
  cat.cex=0.70, main.pos=c(0.5,1.10), cex=0.65, main.cex=0.9
)

venn2 <- venn.diagram(x=list(
  "Breast" = breast.hypo$Probes,
  "Cerebellum" = cerebellum.hypo$Probes,
  "Frontal Cortex" = frontal.hypo$Probes), file=NULL,
  scaled = FALSE,fill = viridis(3), main="Hypovariable Gene Sets", cat.dist=0.09,
  cat.cex=0.70, main.pos=c(0.5,1.10), cex=0.65, main.cex=0.9
)

venn3 <- venn.diagram(x=list(
  "Breast" = breast.NV$Probes,
  "Cerebellum" = cerebellum.NV$Probes,
  "Frontal Cortex" = frontal.NV$Probes), file=NULL,
  scaled = FALSE,fill = viridis(3), main="Non-Variable Gene Sets", cat.dist=0.09,
  cat.cex=0.70, main.pos=c(0.5,1.10), cex=0.65, main.cex=0.9
)



# Chromosomal Map

df <- read.table("Data/Outputs/Final All 3.txt")
breast <- read.table("Data/OUtputs/Breast P1.txt")
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
df$Class[df$Sig == "NS"] <- "NV"
#df <- subset(df, df$Sig == "Sig")


#df <- merge(df, unique(genelist), by=c("Probes", "Gene.Symbol"))
#df$Tissue <- as.factor(df$Tissue)
#df$Class <- as.factor(df$Class)

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
  #scale_color_grey(start=0.5, end=0) +
  #scale_color_viridis(discrete=TRUE, option='A', begin=0.1, end=0.6, direction = -1) +
  scale_color_manual(values= c("#bc2323", '#057010')) +
  guides(color=FALSE)

p3

layout.mat <- matrix(c(1,rep(3,20), 2,rep(4,20), 2,rep(4,20)),
                     ncol=21,byrow=T)

g <- arrangeGrob(arrangeGrob(rectGrob(gp=gpar(col=NA)), top="A"), arrangeGrob(rectGrob(gp=gpar(col=NA)), top="B"),
                 arrangeGrob(gTree(children=venn1), gTree(children=venn2),gTree(children=venn3), ncol=3), 
                 arrangeGrob(p3, vp=viewport(width=1, height=.95)), 
                 layout_matrix=layout.mat, vp=viewport(width=0.95, height=.95))
grid.draw(g)

ggsave("Figures/Figure 4.tiff", plot=g, width=7.5, height=(8.75/1.5), units = 'in')

