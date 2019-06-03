setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
library(dplyr)

ggDF <- read.table("Data/Outputs/Final All 3.txt")
ggDF$Class <- revalue(ggDF$Class, c("NV" = "Non-Variable", "Hypervariable" = "Hyper-Variable",
                                    "Hypovariable" = "Hypo-Variable"))

p1 <- ggplot(ggDF, aes(x=Class, y=Means, fill=Class)) +
  geom_boxplot(width=0.4, outlier.alpha = 0) + facet_grid(.~Tissue) +
  scale_fill_grey(start=0.95, end=0.5) + guides(fill=FALSE) +
  theme(axis.text.y=element_text(size=8),
        axis.title.x = element_text(vjust=5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = -25, vjust=0.5)) +
  ylab("Probability of Matching\nEV Classifications")

p1


#####
# B #
#####


dat <- read.table("Data/Outputs/Final All 3.txt")
dat$Class <- revalue(dat$Class, c("NV" = "Non-Variable", "Hypervariable" = "Hyper-Variable",
                                  "Hypovariable" = "Hypo-Variable"))


sig <- as.data.frame(table(dat$Class, dat$Sig, dat$Tissue))
sig <- subset(sig, sig$Var2 == "Sig")
colnames(sig) <- c("EV", "Geneset", "Tissue", "Count")
sig$Geneset <- "After Retest"

totals <- as.data.frame(table(dat$Class, dat$Tissue))
colnames(totals) <- c("EV", "Tissue", "Count")
totals$Geneset <- "Before Retest"
totals <- totals[,c(1,4,2,3)]

ggDF <- as.data.frame(rbind(sig, totals))
ggDF$Geneset <- factor(ggDF$Geneset, levels=c("Before Retest", "After Retest"))
ggDF <- ggDF[-c(3,6,9),]


# Calculate NV
dat[which(dat$Sig != "Sig"),"Class"] <- "Non-Variable"
NV_df <- dat %>% filter(Class == "Non-Variable")
sig2 <- as.data.frame(table(NV_df$Class, NV_df$Tissue))
sig2$Geneset <- "After Retest"
sig2 <- sig2[,c(1,4,2,3)]
colnames(sig2) <- c("EV", "Geneset", "Tissue", "Count")
sig2 <- sig2 %>% filter(Count > 0)


ggDF <- rbind(ggDF, sig2)


p2 <- ggplot(ggDF, aes(x=Tissue, y=Count, fill=Geneset)) +
  geom_bar(stat='identity', position=position_dodge(width=0.5), width=0.5) + 
  facet_wrap(~EV, scales='free') +
  theme(axis.text.y=element_text(size=8),
        axis.title.x = element_text(vjust=5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = -25, vjust=0.75, hjust=0.2)) +
  scale_fill_grey(start=0.75, end=0.35)
p2



layout.mat <- matrix(c(1,rep(3,20), 1,rep(3,20), 2,rep(4,20), 2,rep(4,20)),
                     ncol=21,byrow=T)

g <- arrangeGrob(arrangeGrob(rectGrob(gp=gpar(col=NA)), top="A"), arrangeGrob(rectGrob(gp=gpar(col=NA)), top="B"),
                 arrangeGrob(p1, vp=viewport(width=1, height=.95)), 
                 arrangeGrob(p2, vp=viewport(width=1, height=.95)), 
                 layout_matrix=layout.mat, vp=viewport(width=0.95, height=.95))
grid.draw(g)

ggsave("Figures/Figure 3.tiff", plot=g, width=7.5, height=(8.75/1.5), units = 'in')
