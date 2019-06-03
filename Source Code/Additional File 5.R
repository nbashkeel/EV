setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(gplots)
pdf("Figures/Additional File 5.pdf", width=16)
textplot("Additional File 5. Complete list of GO terms for essential genes
         ", halign="center", valign='top')

source("Data/Outputs/GO/Essentiality/Treemap/Hypervariable Essential Common.r")
source("Data/Outputs/GO/Essentiality/Treemap/Hypervariable Essential Breast.r")
source("Data/Outputs/GO/Essentiality/Treemap/Hypervariable Essential Cerebellum.r")
source("Data/Outputs/GO/Essentiality/Treemap/Hypervariable Essential Frontal.r")

source("Data/Outputs/GO/Essentiality/Treemap/Hypovariable Essential Common.r")
source("Data/Outputs/GO/Essentiality/Treemap/Hypovariable Essential Breast.r")
source("Data/Outputs/GO/Essentiality/Treemap/Hypovariable Essential Cerebellum.r")
source("Data/Outputs/GO/Essentiality/Treemap/Hypovariable Essential Frontal.r")
dev.off()