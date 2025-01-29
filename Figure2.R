################################################################################
# Author: Addison Dale Buxton-Martin
# 
# Description: This script is associated with the publication ___
# This script will produce the figures for this publication. This script is
# reliant on data already being loaded and helper functions being defined.
#
################################################################################

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


options(contrasts=c("contr.sum", "contr.poly"))

#Running prerequisite scripts if they have not already been run
if(!exists("my_data")){source("LoadData.R")}
if(!exists("FE_analysis")){source("AnalysisFunctions.R")}
if(!exists("Gall_FE")){source("Analysis.R")}

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
  "ggplot2", "tidyverse", "lme4", "lmerTest", "glmmTMB", "DHARMa", "car", "corrplot", "emmeans", "ggpp", "gridExtra", "scales"))


block_formula <- "+ Block.Position1 + Block.Shelf + Block.RowSplit"

emm_vc_nodules <- emmeans(specs = c("RGN"), type = "response",
                                    FE_analysis(my_data %>% mutate("RGN" = interaction(Rhizo, Genotype, Nema)),
                                                formula = "Total.Nodules ~ Volume.mm3 + RGN",
                                                singular.ok = TRUE, plot = FALSE,
                                                return_model = TRUE)) %>%
  as.data.frame() %>% rowwise() %>%
  mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RGN)))[1],
         Genotype = unlist(strsplit(split = "\\.", as.character(RGN)))[2],
         Nema = unlist(strsplit(split = "\\.", as.character(RGN)))[3])

emm_vc_nodules_by_rhizo <- emmeans(specs = c("RN"), type = "response",
                                    FE_analysis(my_data %>% mutate("RN" = interaction(Rhizo, Nema)),
                                                formula = "Total.Nodules ~ Volume.mm3 + RN",
                                                singular.ok = TRUE, plot = FALSE, #block = FALSE,
                                                return_model = TRUE)) %>%
  as.data.frame() %>% rowwise() %>%
  mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RN)))[1],
         Nema = unlist(strsplit(split = "\\.", as.character(RN)))[2])

emm_vc_nodules_by_geno <- emmeans(specs = c("GN"), type = "response",
                                             FE_analysis(my_data %>% mutate("GN" = interaction(Genotype, Nema)),
                                                         formula = "Total.Nodules ~ Volume.mm3 + GN",
                                                         singular.ok = TRUE, plot = FALSE, #block = FALSE,
                                                         return_model = TRUE)) %>%
  as.data.frame() %>% rowwise() %>%
  mutate(Genotype = unlist(strsplit(split = "\\.", as.character(GN)))[1],
         Nema = unlist(strsplit(split = "\\.", as.character(GN)))[2])




p_nodules_geno <- ggplot(emm_vc_nodules_by_geno) + 
  aes(x = Nema, y = response, ymax = response + SE, ymin = response - SE) +
  aes(fill = Genotype) + geom_pointrange(position = position_dodge(width = 0.33)) + 
  geom_line(aes(group = Genotype), position = position_dodge(width = 0.33)) + 
  scale_x_discrete(labels = c("-" = "Uninfected", "+" = "Infected")) +
  xlab("Nematode infection status") +
  ylab("Gall counts\n(volume corrected)") +
  ylim(c(0,110)) + theme_minimal() +
  ggtitle("A") + theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot")

p_nodules_rhizo <- ggplot(emm_vc_nodules_by_rhizo) + 
  aes(x = Nema, y = response, ymax = response + SE, ymin = response - SE) +
  aes(fill = Rhizo) + geom_pointrange(shape = 21, position = position_dodge(width = 0.33)) + 
  geom_line(aes(group = Rhizo), position = position_dodge(width = 0.33)) +
  scale_x_discrete(labels = c("-" = "Uninfected", "+" = "Infected")) +
  xlab("Nematode infection status") +
  ylab("Gall counts\n(volume corrected)") +
  ylim(c(0,110)) + theme_minimal() +
  ggtitle("B") + theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot")


p_nods_variance <- ggplot(Nods_vc_RE %>% filter(`Chi Df` == 1)%>%
                            filter(Effect != "(Intercept)" & Effect != "Block.Position1" & Effect != "scale(Volume.mm3)")) +
  aes(fill = Effect, y = Variance, x = `Chi Df`, label = sig) +
  geom_bar(position = "fill", stat="identity") +
  geom_text(data = . %>%
              filter(sig != "" & sig != ".") %>% 
              mutate(pct = Variance / sum(Variance, na.rm = TRUE)),
            aes(label = paste(percent(pct), sig)),
            position = position_fill(vjust =0.5)) +
  theme_minimal() +
  ggtitle("C") + theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())

p_nods <- grid.arrange(p_nodules_geno, p_nodules_rhizo, p_nods_variance, 
                       ncol = 3, widths = c(3,3,2))

