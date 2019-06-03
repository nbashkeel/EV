setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(rowr)
library(plyr)
library(mixtools)
library(reshape2)
library(mixtools)
library(illuminaio)
library(limma)
library(ggplot2)
library(grid)
library(gridExtra)

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
cere <- read.table("Data/Outputs/Cerebellum Bimodal.txt")
uni <- read.table("Data/Outputs/Cerebellum Unimodal.txt")

# Bimodal
set.seed(1234)
dat <- cerebellum[which(rownames(cerebellum)=="ILMN_1783142"),]
dens <- na.omit(unlist(dat))
dens <- dens[order(dens)]
two.dist <- normalmixEM(dens,k=2,mu=c(min(dens),max(dens)),sigma=c(sd(dens),sd(dens)),arbvar=F,arbmean=F)
EM <- two.dist
x       <- with(EM,seq(min(x),max(x),len=1000))
pars    <- with(EM,data.frame(comp=colnames(posterior), mu, sigma,lambda))
em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))
em.df[nrow(em.df) + 1,] = list(6.426527,'comp.1',0,0,0,0)
em.df[nrow(em.df) + 1,] = list(14.08107,'comp.2',0,0,0,0)

bimodal <- ggplot(data.frame(x=EM$x),aes(x,y=..density..)) + 
  geom_line(stat="density", linetype=2, size=0.75) +
  geom_point(data=em.df[1:(nrow(em.df)-2),],aes(x,y,color=comp),size=0.5, show.legend = F)+
  geom_polygon(data=em.df,aes(x,y,fill=comp),alpha=0.25, show.legend = F)+
  scale_color_manual(values= c('#fc8d59', "#91bfdb")) +
  theme_bw() +
  xlab("Expression") + ylab("Density") + ggtitle("Bimodal Expression Distribution") +
  geom_vline(xintercept = pars$mu, linetype=3) +
  theme(axis.text.y=element_text(size=8),
        axis.title.x = element_text(vjust=1), plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))


# Unimodal
set.seed(1234)
dat <- cerebellum[which(rownames(cerebellum)=="ILMN_2162860"),]
dens <- na.omit(unlist(dat))
dens <- dens[order(dens)]
two.dist <- normalmixEM(dens,k=2,mu=c(min(dens),max(dens)),sigma=c(sd(dens),sd(dens)),arbvar=F,arbmean=F)
EM <- two.dist
x       <- with(EM,seq(min(x),max(x),len=1000))
pars    <- with(EM,data.frame(comp=colnames(posterior), mu, sigma,lambda))
em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))

em.df[nrow(em.df) + 1,] = list(7.571703,'comp.1',0,0,0,0)
em.df[nrow(em.df) + 1,] = list(9.122956,'comp.2',0,0,0,0)



unimodal <- ggplot(data.frame(x=EM$x),aes(x,y=..density..)) + 
  geom_line(stat="density", linetype=2, size=0.75) +
  geom_point(data=em.df[1:(nrow(em.df)-4),],aes(x,y,color=comp),size=0.5, show.legend = F)+
  geom_polygon(data=em.df,aes(x,y,fill=comp),alpha=0.25, show.legend = F)+
  scale_color_manual(values= c('#fc8d59', "#91bfdb")) +
  theme_bw() +
  xlab("Expression") + ylab("Density") + ggtitle("Unimodal Expression Distribution") +
  geom_vline(xintercept = pars$mu, linetype=3) +
  theme(axis.text.y=element_text(size=8),
        axis.title.x = element_text(vjust=1), plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))

p3 <- grid.arrange(bimodal, unimodal, ncol=2)
p3


ggsave("Figures/Figure 2.tiff", plot=p3, width=7.5, height=8.75/3, units="in")
