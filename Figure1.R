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
if(!exists("my_data")){source(file.path(dir_main, "LoadData.R"))}
if(!exists("FE_analysis")){source(file.path(dir_main, "AnalysisFunctions.R"))}
if(!exists("Gall_FE")){source(file.path(dir_main, "Analyses.R"))}

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
  "ggplot2", "tidyverse", "lme4", "lmerTest", "glmmTMB", "DHARMa", "car", "corrplot", "emmeans", "ggpp", "gridExtra", "ggpubr", "scales"))


block_formula <- "+ Block.Position1 + Block.Shelf + Block.RowSplit"


emm_vc_galls <- emmeans(specs = c("RG"), type = "response",
                                  FE_analysis(my_data %>% filter(Nema == "+") %>% mutate("RG" = interaction(Rhizo,Genotype)),
                                              formula = "Galls ~ Volume.mm3 + RG",
                                              singular.ok = TRUE, plot = FALSE, #block = FALSE,
                                              return_model = TRUE)) %>%
  as.data.frame() %>% rowwise() %>%
  mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RG)))[1],
         Genotype = unlist(strsplit(split = "\\.", as.character(RG)))[2])

emm_vc_galls_by_geno <- emmeans(specs = c("Genotype"), type = "response",
                                  FE_analysis(my_data %>% filter(Nema == "+"),
                                              formula = "Galls ~ Volume.mm3 + Genotype",
                                              singular.ok = TRUE, plot = FALSE, #block = FALSE,
                                              return_model = TRUE)) %>% as.data.frame()

emm_vc_galls_by_rhizo <- emmeans(specs = c("Rhizo"), type = "response",
                                          FE_analysis(my_data %>% filter(Nema == "+"),
                                                      formula = "Galls ~ Volume.mm3 + Rhizo",
                                                      singular.ok = TRUE, plot = FALSE, #block = FALSE,
                                                      return_model = TRUE)) %>% as.data.frame()


emm_vc_galls_subclusters <- rbind(
  emmeans(specs = c("RG"), type = "response",
          FE_analysis(subcluster_1 %>% filter(Nema == "+") %>% mutate("RG" = interaction(Rhizo,Genotype)),
                      formula = "Galls ~ Volume.mm3 + RG",
                      singular.ok = TRUE, block = FALSE, plot = FALSE,
                      return_model = TRUE)) %>%
    as.data.frame() %>% rowwise() %>%
    mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RG)))[1],
           Genotype = unlist(strsplit(split = "\\.", as.character(RG)))[2],
           subcluster = "Subcluster 1"),
  emmeans(specs = c("RG"), type = "response",
          FE_analysis(subcluster_2 %>% filter(Nema == "+") %>% mutate("RG" = interaction(Rhizo,Genotype)),
                      formula = "Galls ~ Volume.mm3 + RG",
                      singular.ok = TRUE, block = FALSE, plot = FALSE,
                      return_model = TRUE)) %>%
    as.data.frame() %>% rowwise() %>%
    mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RG)))[1],
           Genotype = unlist(strsplit(split = "\\.", as.character(RG)))[2],
           subcluster = "Subcluster 2"),
  emmeans(specs = c("RG"), type = "response",
          FE_analysis(subcluster_3 %>% filter(Nema == "+") %>% mutate("RG" = interaction(Rhizo,Genotype)),
                      formula = "Galls ~ Volume.mm3 + RG",
                      singular.ok = TRUE, block = FALSE, plot = FALSE,
                      return_model = TRUE)) %>%
    as.data.frame() %>% rowwise() %>%
    mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RG)))[1],
           Genotype = unlist(strsplit(split = "\\.", as.character(RG)))[2],
           subcluster = "Subcluster 3"),
  emmeans(specs = c("RG"), type = "response",
          FE_analysis(subcluster_4 %>% filter(Nema == "+") %>% mutate("RG" = interaction(Rhizo,Genotype)),
                      formula = "Galls ~ Volume.mm3 + RG",
                      singular.ok = TRUE, block = FALSE, plot = FALSE,
                      return_model = TRUE)) %>%
    as.data.frame() %>% rowwise() %>%
    mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RG)))[1],
           Genotype = unlist(strsplit(split = "\\.", as.character(RG)))[2],
           subcluster = "Subcluster 4"),
  emmeans(specs = c("RG"), type = "response",
          FE_analysis(subcluster_5 %>% filter(Nema == "+") %>% mutate("RG" = interaction(Rhizo,Genotype)),
                      formula = "Galls ~ Volume.mm3 + RG",
                      singular.ok = TRUE, block = FALSE, plot = FALSE,
                      return_model = TRUE)) %>%
    as.data.frame() %>% rowwise() %>%
    mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RG)))[1],
           Genotype = unlist(strsplit(split = "\\.", as.character(RG)))[2],
           subcluster = "Subcluster 5"))





default_aes <- aes(y = response, ymin = response - SE, ymax = response + SE)

p_galls_by_geno <- ggplot(emm_vc_galls_by_geno) +
  aes(x = Genotype) + default_aes + 
  geom_pointrange(size = 0.25) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ylab("Gall counts\n(volume corrected)") +
  xlab("Grouped by plant genotype") +
  ylim(c(0,70)) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5))

p_galls_by_rhizo <- ggplot(emm_vc_galls_by_rhizo) +
  aes(x = Rhizo) + default_aes +
  geom_pointrange(size = 0.25) +
  geom_hline(yintercept = 0) + 
  theme_minimal() +
  xlab("Grouped by Rhizobia strain") +
  scale_y_continuous(position = "right", limits = c(0,70)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))

p_galls_genorhizo <- ggarrange(p_galls_by_geno, p_galls_by_rhizo, ncol = 2, align = c("h")) +
  ggtitle("A") + theme(plot.title = element_text(size = 18, hjust = 0)) +
  theme(plot.title.position = "plot")

p_galls_by_genorhizo <- ggplot(emm_vc_galls) +
  default_aes + aes(x = Genotype, fill = Rhizo) + 
  geom_pointrange(shape = 21, position = position_dodge(width = 0.7)) +
  geom_vline(xintercept = seq(1.5,16.5,1), color = "darkgray") +
  geom_hline(yintercept = c(-0.1,0)) +
  theme_minimal() + theme(panel.grid.major.x = element_blank()) +
  ylab("Gall counts\n(volume corrected)") +
  ylim(c(0,70)) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
  theme(plot.title.position = "plot") +
  ggtitle("B") + theme(plot.title = element_text(size = 18, hjust = 0))

p_galls_subclusters <- ggplot(emm_vc_galls_subclusters) +
  aes(x = Genotype, fill = Rhizo) + default_aes +
  geom_line(aes(group = Rhizo, color = Rhizo), position = position_dodge(width = 0.25)) + 
  geom_pointrange(shape = 21, position = position_dodge(width = 0.4)) +
  theme_classic() +
  geom_text(data = as.data.frame(rbind(c("HM049", 60, "", "Subcluster 1"), c("HM012", 60, "", "Subcluster 2"), c("HM011", 60, "", "Subcluster 3"),
                                       c("HM012", 60, "*", "Subcluster 4"), c("HM071", 60, "*", "Subcluster 5"))) %>%
              mutate(V2 = as.numeric(V2), subcluster = V4),
            aes(x = V1, y = V2, label = V3), inherit.aes = FALSE,
            position = position_nudge(x = 0.5), size = 5) +
  facet_wrap(~ subcluster, scales = "free_x", nrow = 1) +
  ylab("Gall counts\n(volume corrected)") +
  ylim(c(0,70)) +
  theme(plot.title.position = "plot") +
  ggtitle("C") + theme(plot.title = element_text(size = 18, hjust = 0))
p_galls_subclusters

p_galls_variance <- ggplot(
  Gall_vc_RE %>%
    filter(`Chi Df` == 1) %>%
    filter(Effect != "(Intercept)" & Effect != "Block.Position1" & Effect != "scale(Volume.mm3)")) +
  aes(fill = Effect, y = Variance, x = `Chi Df`) +
  geom_bar(position = "fill", stat="identity") +
  geom_text(data = . %>%
              filter(sig != "" & sig != ".") %>% 
              mutate(pct = Variance / sum(Variance, na.rm = TRUE)),
            aes(label = paste(percent(pct), sig)),
            position = position_fill(vjust =0.5)) +
  theme_minimal() + theme(legend.position = "bottom") +
  ggtitle("D") + theme(plot.title = element_text(size = 18, hjust = 0)) +
  theme(plot.title.position = "plot") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  guides(fill=guide_legend(ncol = 1))

p_galls <- grid.arrange(p_galls_genorhizo, p_galls_by_genorhizo, p_galls_subclusters, p_galls_variance,
                     layout_matrix = rbind(c(1,4),c(2,4),c(3,4)),
                     widths = c(4,1))
