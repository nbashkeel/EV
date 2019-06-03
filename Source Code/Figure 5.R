setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(mixtools)
library(ggplot2)
library(grid)
library(gridExtra)
library(gridGraphics)
library(reshape2)
library(dplyr)

dat <- read.table("Data/Outputs/Final Meth Summary.txt")


ggDF2 <- subset(dat, dat$Class != "NV")


########################
# Methylation Clusters #
########################

cere.dens <- subset(dat, dat$Tissue == "Cerebellum")
cere.meth.range <- range(dat$Meth.Median[dat$Methylation.Cluster == "Medium Methylation"])
cere.dens <- density(cere.dens$Meth.Median, n=nrow(dat), from=0,to=1)
cere.dens <- data.frame(x=cere.dens$x, y=cere.dens$y)
cere.dens <- cere.dens[order(cere.dens$x),]
cere.dens$Methylation.Cluster <- "Medium Methylated"
cere.dens$Methylation.Cluster[cere.dens$x < cere.meth.range[1]] <- "Non-Methylated"
cere.dens$Methylation.Cluster[cere.dens$x > cere.meth.range[2]] <- "Highly Methylated"
cere.dens$Tissue <- "Cerebellum"

frontal.dens <- subset(dat, dat$Tissue == "Frontal Cortex")
frontal.meth.range <- range(dat$Meth.Median[dat$Methylation.Cluster == "Medium Methylation"])
frontal.dens <- density(frontal.dens$Meth.Median, n=nrow(dat), from=0,to=1)
frontal.dens <- data.frame(x=frontal.dens$x, y=frontal.dens$y)
frontal.dens <- frontal.dens[order(frontal.dens$x),]
frontal.dens$Methylation.Cluster <- "Medium Methylated"
frontal.dens$Methylation.Cluster[frontal.dens$x < frontal.meth.range[1]] <- "Non-Methylated"
frontal.dens$Methylation.Cluster[frontal.dens$x > frontal.meth.range[2]] <- "Highly Methylated"
frontal.dens$Tissue <- "Frontal Cortex"

ggDF <- rbind(cere.dens, frontal.dens)
ggDF$Methylation.Cluster <- factor(ggDF$Methylation.Cluster, levels=c("Non-Methylated", "Medium Methylated", "Highly Methylated"))
ggDF$Class <- revalue(ggDF$Class, c("NV" = "Non-Variable", "Hypervariable" = "Hyper-Variable",
                                    "Hypovariable" = "Hypo-Variable"))


p1 <- ggplot(ggDF, aes(x,y)) + 
  geom_ribbon(aes(ymin=0, ymax=y, fill=Methylation.Cluster)) + 
  geom_line(alpha=0.5) +
  scale_fill_brewer(name="Methylation Cluster", type='div', palette=7) +
  xlab("Methylation") + ylab("\nDensity") +
  facet_grid(.~Tissue) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.spacing = unit(1, "lines"),
        axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0), breaks=c(seq(0,1,0.25)), labels=c(0,0.25,0.50,0.75,1.0))

p1




dummy2 <- by(ggDF2, ggDF2$Tissue, function(x){
  range(subset(x, x$Methylation.Cluster == "Medium Methylation")$Meth.Median)
})
dummy2 <- do.call(rbind,dummy2)
dummy2 <- melt(dummy2, id.vars="row.names")
colnames(dummy2) <- c("Tissue", "x", "Cutoff")

dat$Class <- revalue(dat$Class, c("NV" = "Non-Variable", "Hypervariable" = "Hyper-Variable",
                                  "Hypovariable" = "Hypo-Variable"))

p2 <- ggplot(dat, aes(x=Meth.Median, color = Class)) +
  geom_line(stat='Density') + facet_grid(.~Tissue) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_vline(data = dummy2, aes(xintercept = Cutoff), alpha=0.55, lty=5) +
  xlab("Median Methylation") + ylab("Density") +
  theme(panel.spacing = unit(1, "lines"),
        axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0), breaks=c(seq(0,1,0.25)), labels=c(0,0.25,0.50,0.75,1.0)) +
  scale_color_manual(values=c("#2F6C8EFF", "#bc2323", "#38B977FF"))


p2


layout.mat <- matrix(c(1,rep(3,20), 2,rep(4,20)),
                     ncol=21,byrow=T)

g <- arrangeGrob(arrangeGrob(rectGrob(gp=gpar(col=NA)), top="A"), arrangeGrob(rectGrob(gp=gpar(col=NA)), top="B"), 
                 arrangeGrob(p1,vp=viewport(width=1, height=.95)), 
                 arrangeGrob(p2,vp=viewport(width=1, height=.95)), 
                 layout_matrix=layout.mat, vp=viewport(width=0.95, height=.95))
grid.draw(g)
ggsave("Figures/Figure 5.tiff", plot=g, width=7.5, height=(8.75/2), units="in")
