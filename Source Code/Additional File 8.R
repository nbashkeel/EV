setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(plyr)
library(ggplot2)
library(gplots)
library(grid)
library(gridExtra)

df <- read.table("C:/Users/Bashk/Desktop/Thesis\ Backup/2018/Data/Full Exp Meth DF.txt")
df$MV.Class <- revalue(df$MV.Class, c("Hypervariable" = "Hyper-MV",
                                          "Hypovariable" = "Hypo-MV",
                                          "NV" = "Non-MV"))
df$EV.Class <- revalue(df$EV.Class, c("NV" = "Non-Variable"))


cere <- subset(df, df$Tissue == "Cerebellum")
frontal <- subset(df, df$Tissue == "Frontal Cortex")



p1 <- ggplot(cere, aes(x=EV, y=MV)) +
  geom_point(cex=0.5) + facet_wrap(MV.Class~EV.Class, scales='free')
p2 <- ggplot(frontal, aes(x=EV, y=MV)) +
  geom_point(cex=0.5) + facet_wrap(MV.Class~EV.Class, scales='free')


EV_MV_cere <- ddply(cere, .(EV.Class, MV.Class), function(x) cor(x$EV, x$MV, method='spearman'))
colnames(EV_MV_cere) <- c("EV Class", "MV Class", "Spearman Correlation")
EV_MV_fc <- ddply(frontal, .(EV.Class, MV.Class), function(x) cor(x$EV, x$MV, method='spearman'))
colnames(EV_MV_fc) <- c("EV Class", "MV Class", "Spearman Correlation")



pdf("Figures/Additional File 8.pdf")
textplot("Additional File 8. Methylation Variability
", halign="center", valign='top')
grob(textplot(EV_MV_cere, show.rownames = FALSE))
title("Correlation between EV and MV in Cerebellum Tissue")
p1

grob(textplot(EV_MV_fc, show.rownames = FALSE))
title("Correlation between EV and MV in Frontal Cortex Tissue")
p2
dev.off()
