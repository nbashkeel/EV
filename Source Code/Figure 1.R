setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(VennDiagram)
library(plyr)
library(base2grob)
library(plotrix)

color_palette <- c('#fc8d59', "#ffffbf", "#91bfdb")

breast <- read.table("Data/Outputs/Breast P1.txt")
breast$Tissue <- "Breast"
cerebellum <- read.table("Data/Outputs/Cerebellum P1.txt")
cerebellum$Tissue <- "Cerebellum"
frontal <- read.table("Data/Outputs/Frontal P1.txt")
frontal$Tissue <- "Frontal Cortex"

df <- rbind(breast, cerebellum, frontal)
df$Tissue <- as.factor(df$Tissue)

p1 <- ggplot(df, aes(x=Median, y=Obs.MAD)) +
  stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200)+
  scale_alpha_manual(values=c(0,1),guide="none") + 
  facet_wrap(~Tissue, scales="free_x") + ylim(0,1.2) +
  scale_fill_gradient2(low="#ffffff", high="#000000", midpoint=0.75) +
  xlab("Median Expression") + ylab("MAD") +
  guides(fill=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_point(alpha=0.1, size=0.01) +
  geom_line(aes(x=Median, y=Exp.MAD),col='red',alpha=0.75)
p1


EV.df <- c(breast$EV, cerebellum$EV,frontal$EV)
EV.df <- data.frame(EV=EV.df, Class=c(rep("Breast", nrow(breast)),
                                      rep("Cerebellum",nrow(cerebellum)),
                                      rep("Frontal", nrow(frontal))))

abline.df <- data.frame(Class=c("Breast", "Cerebellum", "Frontal"),
                        llim=c(min(breast$EV[breast$Class=="NV"]),
                               min(cerebellum$EV[cerebellum$Class=="NV"]),
                               min(frontal$EV[frontal$Class=="NV"])),
                        ulim=c(max(breast$EV[breast$Class=="NV"]),
                               max(cerebellum$EV[cerebellum$Class=="NV"]),
                               max(frontal$EV[frontal$Class=="NV"])))

abline.df <- melt(abline.df, id.vars = "Class")

f2a_func <- function(){
  d <- density(breast$EV)
  par(mar=c(4.0,3.0,3,0))
  par(las=1)
  par(mgp=c(2,0.5,0))
  par(cex.axis=0.75, cex.main=0.75)
  par(tck=-0.01)
  gap.plot(c(-0.5,d$x), c(0,d$y), gap=range(c(0.7,2.1)), 
           gap.axis="x", type="l", xlab="", ylim=c(0,12),
           ylab="Density", 
           xtics=c(seq(-0.5,0.5,by=0.5),2.25),
           ytics=c(seq(0,ceiling(max(d$y)), length.out=4)))
  title("Breast", line=0.25, font.main = 1, cex=1)
  polygon(d, col=color_palette[1], border="black")
  abline(v=abline.df[c(1,4),3])}

f2a <- base2grob(f2a_func)

f2b_func <- function(){
  d <- density(cerebellum$EV)
  par(mar=c(4.0,3.0,3,0))
  par(las=1)
  par(mgp=c(2,0.5,0))
  par(cex.axis=0.75, cex.main=0.75)
  par(tck=-0.01)
  gap.plot(c(-0.2,d$x), c(0,d$y), gap=range(c(0.29,0.8)), 
           gap.axis="x", type="l", xlab="Expression Variability", 
           ylab="", ylim=c(0,18.25),
           xtics=c(seq(-0.2,0.4,by=0.2),0.85),
           ytics=c(seq(0,ceiling(max(d$y)), length.out=4)))
  title("Cerebellum", line=0.25, font.main = 1, cex=1)
  polygon(d, col=color_palette[2], border="black")
  abline(v=abline.df[c(2,5),3])}

f2b <- base2grob(f2b_func)

f2c_func <- function(){
  d <- density(frontal$EV)
  par(mar=c(4.0,3.0,3,0))
  par(las=1)
  par(mgp=c(2,0.25,0))
  par(mgp=c(1,0.5,0))
  par(cex.axis=0.75, cex.main=0.75)
  par(tck=-0.01)
  gap.plot(c(-0.2,d$x), c(0,d$y), gap=range(c(0.29,1.1)), 
           gap.axis="x", type="l", xlab="", ylim=c(0,18),
           ylab="",
           xtics=c(seq(-0.2,0.4,by=0.2),1.15),
           ytics=c(0,6,12,18))
  title("Frontal Cortex", line=0.25, font.main = 1, cex=1)
  polygon(d, col=color_palette[3], border="black")
  abline(v=abline.df[c(2,5),3])}

f2c <- base2grob(f2c_func)

p2 <- arrangeGrob(f2a, f2b, f2c, ncol=3)


df <- rbind(breast, cerebellum, frontal)
df$Tissue <- as.factor(df$Tissue)

ggplotRegression <- function (dat) {
  fit <- lm(EV ~ Median, data = dat)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(size=0.25, alpha=0.5) +
    stat_smooth(method = "lm", col = "red", lwd=0.5) +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "\nIntercept =",signif(fit$coef[[1]],5 ),
                       "\nSlope =",signif(fit$coef[[2]], 5),
                       "\nP=",signif(summary(fit)$coef[2,4], 5)),
         cex=0.015)
}

p3a <- ggplotRegression(breast)
p3b <- ggplotRegression(cerebellum)
p3c <- ggplotRegression(frontal)

grid.arrange(p3a, p3b, p3c, ncol = 3)

p3 <- ggplot(df, aes(x=Median, y=EV)) +
  geom_point(alpha=0.25, size=0.5)+
  facet_wrap(~Tissue, scales="free") +
  stat_smooth(method = "lm", col = "red", lwd=0.5) +
  xlab("Median Expression") +
  theme(axis.text.y=element_text(size=8),
        axis.title.x = element_text(vjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))

p3


layout.mat <- matrix(c(1,rep(4,20), 1,rep(4,20),1,rep(4,20),2,rep(5,20),2,rep(5,20),2,rep(5,20), 
                       3,rep(6,20),3,rep(6,20), 3, rep(6,20)),
                     ncol=21,byrow=T)


g <- arrangeGrob(arrangeGrob(rectGrob(gp=gpar(col=NA)), top="A"), arrangeGrob(rectGrob(gp=gpar(col=NA)), top="B"),arrangeGrob(rectGrob(gp=gpar(col=NA)), top="C"),
                 arrangeGrob(p1,vp=viewport(width=1, height=.95)), arrangeGrob(p2,vp=viewport(width=1, height=.95)), arrangeGrob(p3,vp=viewport(width=1, height=.95)), 
                 layout_matrix=layout.mat, vp=viewport(width=0.95, height=.95))
grid.draw(g)

ggsave("Figures/Figure 1.tiff", plot=g, width=7.5, height=8.75, units="in")
