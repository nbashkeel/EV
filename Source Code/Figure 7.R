setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(ggplot2)
library(dplyr)
library(plyr)

methylation <- read.table("Data/Outputs/Final Meth Summary.txt")
ggDF <- methylation %>% filter(Class != 'NV') %>% droplevels()
ggDF$Class <- revalue(ggDF$Class, c("Hypervariable" = "Hyper-Variable",
                                    "Hypovariable" = "Hypo-Variable"))

p3 <- ggplot(ggDF, aes(x=COR)) +
  geom_histogram(binwidth=0.025, col="white") + 
  facet_grid(Class~Tissue) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  xlab("Pearson Correlation") + ylab("Number of Genes") +
  scale_fill_grey(start=0.3, end=0.6) + xlim(-0.4,0.4)
p3

ggsave("Figures/Figure 7.tiff", plot=p3, width=7.5, height=(8.75/2), units="in")
