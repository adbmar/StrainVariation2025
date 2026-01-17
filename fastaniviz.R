################################################################################
# Author: Addison Dale Buxton-Martin
# 
# Description: This script runs the analysis for ___
# This analysis will import data included in its encasing directory
# and perform analyses as outlined by ___
#
################################################################################

options(contrasts=c("contr.sum", "contr.poly"))

###################################
### Loading requisite libraries ###
###################################

# loading requisite libraries
# Define a function to check, install, and load a library
check_and_load_libraries <- function(library_names) {
  # Ensure library_names is a vector
  if (!is.character(library_names)) {
    stop("The input must be a character vector of package names.")
  }
  
  for (library_name in library_names) {
    if (!requireNamespace(library_name, quietly = TRUE)) {
      # Try installing the package using base install.packages
      tryCatch({
        install.packages(library_name, dependencies = TRUE)
      }, error = function(e) {
        # If base install.packages fails, try BiocManager::install
        if (requireNamespace("BiocManager", quietly = TRUE)) {
          BiocManager::install(library_name, dependencies = TRUE)
        } else {
          stop("BiocManager is required but not installed.")
        }
      })
    }
    # Load the library
    library(library_name, character.only = TRUE)
  }
}

check_and_load_libraries(c(
  "ggplot2", "tidyverse", "lme4", "lmerTest", "glmmTMB", "DHARMa", "car", "corrplot", "emmeans", "ggpp"))


################################################################################
############################## DIRECTORY HANDLING ##############################
################################################################################
dir_main <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  dirname(normalizePath(sys.frame(1)$ofile))
}

#manual back up for setting main directory)
################################################################################
# Edit the file path of dir_main to be the folder containing this script.
# The data should be in this main directory as specified in the above Data
# Requirements section in the heading of this script
# dir_phenomain <- file.path("~", "Desktop", "analysis")
################################################################################

#   "GCA_033715635.1_ASM3371563v1_genomic.fna", "MAG177 (D)"
#   "GCA_033716455.1_ASM3371645v1_genomic.fna", "MAG154 (G)"
#   "GCA_033716615.1_ASM3371661v1_genomic.fna", "MAG282 (F)"
#   "GCA_033716655.1_ASM3371665v1_genomic.fna", "MAG358 (I)"
#   "GCA_033714255.1_ASM3371425v1_genomic.fna", "MAG533 (H)"
#   "GCA_033716125.1_ASM3371612v1_genomic.fna", "MAG540 (J)"
#   "GCA_002197105.1_ASM219710v1_genomic.fna", "KH35c (C)"
#   "GCA_009599945.1_ASM959994v1_genomic.fna", "KH48e (E)"
#   "GCA_002197445.1_ASM219744v1_genomic.fna", "USDA1021"
#   "GCA_013315775.1_ASM1331577v1_genomic.fna", "WSM1022"

fastani <- read.csv(file.path(dir_main, "fastanioutput.txt"), sep = "\t", header = FALSE)
colnames(fastani) <- c("query", "reference", "ANI", "bi_frag_maps", "qlength")
fastani <- fastani %>% rowwise() %>%
  mutate(query = gsub("^.*/", "", query),
         reference = gsub("^.*/", "", reference)) %>%
  mutate(query = gsub("GCA_033715635.1_ASM3371563v1_genomic.fna", "MAG177 (D)", query),
         reference = gsub("GCA_033715635.1_ASM3371563v1_genomic.fna", "MAG177 (D)", reference),
         
         query = gsub("GCA_033716455.1_ASM3371645v1_genomic.fna", "MAG154 (G)", query),
         reference = gsub("GCA_033716455.1_ASM3371645v1_genomic.fna", "MAG154 (G)", reference),
         
         query = gsub("GCA_033716615.1_ASM3371661v1_genomic.fna", "MAG282 (F)", query),
         reference = gsub("GCA_033716615.1_ASM3371661v1_genomic.fna", "MAG282 (F)", reference),
         
         query = gsub("GCA_033716655.1_ASM3371665v1_genomic.fna", "MAG358 (I)", query),
         reference = gsub("GCA_033716655.1_ASM3371665v1_genomic.fna", "MAG358 (I)", reference),
         
         query = gsub("GCA_033714255.1_ASM3371425v1_genomic.fna", "MAG533 (H)", query),
         reference = gsub("GCA_033714255.1_ASM3371425v1_genomic.fna", "MAG533 (H)", reference),
         
         query = gsub("GCA_033716125.1_ASM3371612v1_genomic.fna", "MAG540 (J)", query),
         reference = gsub("GCA_033716125.1_ASM3371612v1_genomic.fna", "MAG540 (J)", reference),
         
         query = gsub("GCA_002197105.1_ASM219710v1_genomic.fna", "KH35c (C)", query),
         reference = gsub("GCA_002197105.1_ASM219710v1_genomic.fna", "KH35c (C)", reference),
         
         
         query = gsub("GCA_009599945.1_ASM959994v1_genomic.fna", "KH48e (E)", query),
         reference = gsub("GCA_009599945.1_ASM959994v1_genomic.fna", "KH48e (E)", reference),
         
         
         query = gsub("GCA_002197445.1_ASM219744v1_genomic.fna", "USDA1021", query),
         reference = gsub("GCA_002197445.1_ASM219744v1_genomic.fna", "USDA1021", reference),
         
         query = gsub("GCA_013315775.1_ASM1331577v1_genomic.fna", "WSM1022", query),
         reference = gsub("GCA_013315775.1_ASM1331577v1_genomic.fna", "WSM1022", reference))

similarity_matrix <- reshape2::acast(fastani, query ~ reference, value.var = "ANI", fill = NA)
similarity_matrix[is.na(similarity_matrix)] <- 0
min(similarity_matrix)
distance_matrix <- 1 - similarity_matrix / 100
hc <- hclust(as.dist(distance_matrix), method = "average")
# hc_inv <- hclust(as.dist(similarity_matrix), method = "average")
plot(hc, main = "Dendrogram of FASTANI Results", xlab = "Strain", sub = "", cex = 0.8, ylab = "ANI Distance")
# plot(hc_inv, main = "Dendrogram of FASTANI Results", xlab = "Strain", sub = "", cex = 0.8)
setwd(dir_out)
png(filename = "FASTANIDendrogram.png", width = 2160, height= 2160, res = 400)
plot(hc, main = "Dendrogram of FASTANI Results", xlab = "Strain", sub = "", cex = 0.8, ylab = "ANI Distance")
dev.off()