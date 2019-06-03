setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

pdf("Figures/Additional File 6.pdf", width=16)
textplot("Additional File 6. Complete list of GO terms for age-regulated Hyper-variable genes
         ", halign="center", valign='top')

source("Data/Outputs/GO/Age/REVIGO/Cere Up Gold.r")
source("Data/Outputs/GO/Age/REVIGO/Cere Up Darkorange.r")
source("Data/Outputs/GO/Age/REVIGO/Cere Down Yellow.r")
source("Data/Outputs/GO/Age/REVIGO/Cere Down Green.r")

source("Data/Outputs/GO/Age/REVIGO/Frontal Up Gold.r")
source("Data/Outputs/GO/Age/REVIGO/Frontal Up Darkorange.r")
source("Data/Outputs/GO/Age/REVIGO/Frontal Up Yellow.r")
source("Data/Outputs/GO/Age/REVIGO/Frontal Down Red.r")

source("Data/Outputs/GO/Age/REVIGO/Cere Up Darkorange Location.r")
source("Data/Outputs/GO/Age/REVIGO/Cere Down Yellow Location.r")
source("Data/Outputs/GO/Age/REVIGO/Cere Down Green Location.r")

source("Data/Outputs/GO/Age/REVIGO/Frontal Up Gold Location.r")
source("Data/Outputs/GO/Age/REVIGO/Frontal Up Darkorange Location.r")
source("Data/Outputs/GO/Age/REVIGO/Frontal Up Yellow Location.r")
source("Data/Outputs/GO/Age/REVIGO/Frontal Down Red Location.r")

dev.off()

