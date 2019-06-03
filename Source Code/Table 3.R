setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(rowr)

df <- read.table("Data/Final All 3.txt")
df$Class[df$Sig == "NS"] <- "NV"
breast <- subset(df, df$Tissue == "Breast")
cerebellum <- subset(df, df$Tissue == "Cerebellum")
frontal <- subset(df, df$Tissue == "Frontal Cortex")

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

Venn.Overlaps <- function(class) {
  breast <- subset(breast, breast$Class == class)
  cerebellum <- subset(cerebellum, cerebellum$Class == class)
  frontal <- subset(frontal, frontal$Class == class)
  
  write.table(data.frame(" "=unique(intersect(intersect(breast$Gene.Symbol, cerebellum$Gene.Symbol),frontal$Gene.Symbol))),
              paste0("Data/Outputs/GO/Overall/Common ", simpleCap(class),".txt"),
              col.names=F, row.names=F, quote = F)
  
  write.table(data.frame(unique(setdiff(setdiff(breast$Gene.Symbol, cerebellum$Gene.Symbol),frontal$Gene.Symbol))),
              paste0("Data/Outputs/GO/Overall/Breast ", simpleCap(class),".txt"),
              col.names=F, row.names=F, quote = F)
  
  write.table(data.frame(unique(setdiff(setdiff(cerebellum$Gene.Symbol, breast$Gene.Symbol),frontal$Gene.Symbol))),
              paste0("Data/Outputs/GO/Overall/Cerebellum ", simpleCap(class),".txt"),
              col.names=F, row.names=F, quote = F)
  
  write.table(data.frame(unique(setdiff(setdiff(frontal$Gene.Symbol, cerebellum$Gene.Symbol),breast$Gene.Symbol))),
              paste0("Data/Outputs/GO/Overall/Frontal ", simpleCap(class),".txt"),
              col.names=F, row.names=F, quote = F)
}

Venn.Overlaps("Hypervariable")
Venn.Overlaps("Hypovariable")
Venn.Overlaps("NV")

# Run gene lists through ConsensusPathDB to get GO Terms: http://cpdb.molgen.mpg.de/
# Run GO Term lists through REVIGO to get treemaps: http://revigo.irb.hr/
  # Save REVIGO Treemaps code to Data/Outputs/GO/Overall/Treemap
  # Run Supplementary Figure 3 Script
