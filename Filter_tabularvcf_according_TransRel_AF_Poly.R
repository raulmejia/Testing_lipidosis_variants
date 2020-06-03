#####################
##### The aim of this program is to apply the following filter to a data frame that contain information about variants 
## AF , Polyphen2 , gene_symbol, TransRel, CHROM, POS
# the rscript way and the manual way doenst retrieve the same result table 41 vs 38 and I dont know why
#####################
# Fix why do line  "Secondary_variant_types_AF_upto_G <- Secondary_variant_types %>% filter(AF <= G) %>% filter(PolyPhen2 > P)" add extra rows 
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
if (!require("limma")) {
  install.packages("limma", dependencies = TRUE)
  library(limma)
}
#############################
### Data given by the user
#############################
#myargs <- commandArgs( trailingOnly = TRUE )
#path_to_your_df <- myargs[1]
#results_path <- myargs[2]
#results_label <- myargs[3]
#G  <- myargs[4]
#P <- myargs[5]

path_to_your_df <-  "/media/rmejia/mountme88/Projects/Phosoholipidosis/Exome_Lipidosis/Results_needed_From_Run1/Testing_variants/Lysoplex/DiExVar_Run1_in_Lysoplex.tsv"
results_path <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/Exome_Lipidosis/Results_needed_From_Run1/Testing_variants/Lysoplex"
results_label <- "_Filtered_by_"
G  <- 0.1
P <- 0.15
############################
### Normalize paths
############################
dir.create( results_path , recursive = T)
normalizePath( results_path )

############################
### Reading the data
############################
your_df <- read.table( path_to_your_df, sep="\t", header = TRUE, stringsAsFactors = FALSE)
dim(your_df)
head(your_df, n=150)
# Splitting in primary and secondary types
Top_variant_type <- your_df %>%
  filter(str_detect(TransRel, 'splice|missense|nonsense|no-start|no-stop|inframe|frameshift'))
dim(Top_variant_type)
head(Top_variant_type , n= 100)

# filtering
Top_variant_type_AF_upto_G <- Top_variant_type %>% filter(AF < G) %>% filter(PolyPhen2 > P)
Top_variant_type_AF_NA <- Top_variant_type[is.na(Top_variant_type[,"AF"]),]
Top_variant_type_Poly_NA <- Top_variant_type[is.na(Top_variant_type[,"PolyPhen2"]),]
Top_variant_type_AF_filtered <- rbind(Top_variant_type_AF_upto_G , Top_variant_type_AF_NA , Top_variant_type_Poly_NA)

# Ordering
Top_variant_type_AF_filtered <- Top_variant_type_AF_filtered[  order(Top_variant_type_AF_filtered[,"POS"]) , ]
Top_variant_type_AF_filtered[,"CHROM"] <- as.numeric(Top_variant_type_AF_filtered[,"CHROM"])
Top_variant_type_AF_filtered <- Top_variant_type_AF_filtered[  order(Top_variant_type_AF_filtered[,"CHROM"]) , ]
str(Top_variant_type_AF_filtered)
dim(Top_variant_type_AF_filtered)

# Saving your results 
path_write_results_Top_priority <- paste0(results_path, "/", removeExt(  basename(  path_to_your_df), sep=".") ,results_label, "_Filtered_by_","Top_priority_","AF_Lesserthan",G,"Polyphen_Biggerthan",P,".tsv" )
write.table( file=path_write_results_Top_priority , Top_variant_type_AF_filtered , sep="\t" )
dim(Top_variant_type_AF_filtered)

Secondary_variant_types <- your_df %>%
  filter(str_detect(TransRel, 'intronic|intergenic|synonymous|UTR'))
dim(Secondary_variant_types)
class(Secondary_variant_types)
head( Secondary_variant_types, n=137)

# filtering
#Secondary_variant_types_AF_upto_G <- Secondary_variant_types %>% filter(AF < G) %>% filter(PolyPhen2 > P)
#Secondary_variant_types_AF_NA <- Secondary_variant_types[is.na(Secondary_variant_types[,"AF"]),]
#Secondary_variant_types_Poly_NA <- Secondary_variant_types[ is.na(Secondary_variant_types[,"PolyPhen2"]) , ]
#Secondary_variant_types_AF_filtered <- rbind(Secondary_variant_types_AF_upto_G , Secondary_variant_types_AF_NA , Secondary_variant_types_Poly_NA)
#dim(Secondary_variant_types_AF_filtered)
# Ordering
#Secondary_variant_types_AF_filtered <- Secondary_variant_types_AF_filtered[  order(Secondary_variant_types_AF_filtered[,"POS"]) , ]
#Secondary_variant_types_AF_filtered[,"CHROM"] <- as.numeric(Secondary_variant_types_AF_filtered[,"CHROM"])
#Secondary_variant_types_AF_filtered <- Secondary_variant_types_AF_filtered[  order(Secondary_variant_types_AF_filtered[,"CHROM"]) , ]
#str(Secondary_variant_types_AF_filtered)
#dim(Secondary_variant_types_AF_filtered)

# Saving your results 
#path_write_results_secondary_priority <- paste0(results_path, "/", removeExt(  basename(  path_to_your_df), sep=".") ,results_label, "_Filtered_by_","Secondary_priority_","AF_Lesserthan",G,"Polyphen_Biggerthan",P,".tsv" )
path_write_results_secondary_priority <- paste0(results_path, "/", removeExt(  basename(  path_to_your_df), sep=".") ,results_label, "_Filtered_by_","Secondary_priority",".tsv" )
write.table( file=path_write_results_secondary_priority , Secondary_variant_types , sep="\t", row.names=FALSE , quote=TRUE)


