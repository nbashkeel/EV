setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

pdf("Figures/Additional File 3.pdf", width=16)
textplot("Additional File 3. Complete list of GO terms for all genes
         ", halign="center", valign='top')

source("Data/Outputs/GO/Overall/Treemap/Hypervariable Common.r")
source("Data/Outputs/GO/Overall/Treemap/Hypovariable Common.r")

source("Data/Outputs/GO/Overall/Treemap/Hypervariable Breast.r")
source("Data/Outputs/GO/Overall/Treemap/Hypervariable Cerebellum.r")
source("Data/Outputs/GO/Overall/Treemap/Hypervariable Frontal.r")

source("Data/Outputs/GO/Overall/Treemap/Hypovariable Breast.r")
source("Data/Outputs/GO/Overall/Treemap/Hypovariable Cerebellum.r")
source("Data/Outputs/GO/Overall/Treemap/Hypovariable Frontal.r")

source("Data/Outputs/GO/Overall/Treemap/Hypervariable Common Location.r")
source("Data/Outputs/GO/Overall/Treemap/Hypovariable Common Location.r")

source("Data/Outputs/GO/Overall/Treemap/Hypervariable Breast Location.r")
source("Data/Outputs/GO/Overall/Treemap/Hypervariable Cerebellum Location.r")
source("Data/Outputs/GO/Overall/Treemap/Hypervariable Frontal Location.r")

source("Data/Outputs/GO/Overall/Treemap/Hypovariable Breast Location.r")
source("Data/Outputs/GO/Overall/Treemap/Hypovariable Cerebellum Location.r")
source("Data/Outputs/GO/Overall/Treemap/Hypovariable Frontal Location.r")
dev.off()