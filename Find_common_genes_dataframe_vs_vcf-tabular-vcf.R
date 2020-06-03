#####################
##### The aim of this program is to indentify the common genes in a VCF or a "tabular-shaped" VCF compared to an external vcf
## Remember to use an input table tab separated, and " separating the fields
## AF , Polyphen2 , gene_symbol, TransRel, CHROM, POS
#####################
#  
#
#
#####################
# Installing & Loading required libraries 
# ###################
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}
if (!require("readr")) {
  install.packages("readr", dependencies = TRUE)
  library(readr)
}
#############################
### Data given by the user
#############################
myargs <- commandArgs( trailingOnly = TRUE )
path_to_your_problemdf <- myargs[1]
path_to_your_tabularvcf <- myargs[2]
results_path <- myargs[3]
results_label <- myargs[4]


############################
### Normalize paths
############################
dir.create( results_path , recursive = T)
normalizePath( results_path )

############################
### Reading the data
############################
your_vcf <- read.table(path_to_your_tabularvcf, sep="\t", header = TRUE, stringsAsFactors = FALSE)
print(path_to_your_tabularvcf)
Problemdf <- read.table( path_to_your_problemdf, sep="\t", skip = 1, header=TRUE )

############################
## matching and saving
############################
positions <- which( your_vcf[ ,"gene_symbol"] %in% Problemdf[ ,"gene_symbol"] )
Genes_in_your_vcf_also_present_in_theproblemdf <- your_vcf[ ,"gene_symbol"][positions]
Variants_in_your_vcf_also_present_in_theproblemdf <- your_vcf[ positions,]
table(Variants_in_your_vcf_also_present_in_theproblemdf[,"gene_symbol"])
unique(Variants_in_your_vcf_also_present_in_theproblemdf[,"gene_symbol"])
table(Variants_in_your_vcf_also_present_in_theproblemdf[,"TransRel"])

path_write_results <- paste0(results_path, "/", results_label,".tsv" )
write.table( file=path_write_results , Variants_in_your_vcf_also_present_in_theproblemdf , sep="\t", row.names=FALSE, quote=TRUE)



