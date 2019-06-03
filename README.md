# Expression Variability
Source code for "Human Gene Expression Variability and its Dependence on Methylation and Aging"


## Install the Dependencies
R Packages
* base2grob
* dplyr
* ggplot2
* gplots
* grid
* gridExtra
* gridGraphics
* gtools
* illuminaio
* limma
* lsa
* mclust
* mixtools
* plotrix
* plyr
* pvclust
* reshape2
* rowr
* seqinr
* tidyr
* VennDiagram
* viridis

<br>

## Download Required Files
> Save all files to Source Code > Data

Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36192
* Series Matrix File (Referred to as `brain_ds`)
  * Subset rows 55-91 (Referred to as `Full Brain Annotations`)
  * Subset rows 55-56, 64-68 (Referred to as `Brain Samples`)

Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36194
* Series Matrix File (Referred to as `Brain Methylation`)
  * Subset rows 46-83 (Referred to as `Methylation Annotations`)

Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8490
* Full table (Referred to as `Methylation Probes`)

Download from: https://www.ebi.ac.uk/ega/datasets/EGAD00010000212
* Dataset (Referred to as `normals_ExpressionMatrix`)

Download from: https://support.illumina.com/downloads/humanht-12_v3_product_files.html
* BGX Format (Referred to as `HumanHt-12_V3_0_R3_11283641_A`)

Download from: https://genome-euro.ucsc.edu/cgi-bin/hgTables
* Download track: Geneid Genes & output format: selected fields from primary and related tables (Referred to as `hgTables Summary`)
* Download track: Geneid Genes & output format: sequence (Referred to as `hgTables Sequences`)

Download from https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003484
* Table s6 (Referred to as `Essential Genes`)

<br>

## External Resources

For GO Term Enrichment Analyses: http://cpdb.molgen.mpg.de/
* Gene set analysis > over-representation analysis

For Treemaps & GO Term Summarization: http://revigo.irb.hr/
* Select a database with GO term sizes: Homo sapiens
