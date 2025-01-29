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
if(!exists("Gall_FE")){source("Analyses.R")}

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
  "ggplot2", "tidyverse", "lme4", "lmerTest", "glmmTMB", "DHARMa", "car", "corrplot", "emmeans", "ggpp", "gridExtra"))


block_formula <- "+ Block.Position1 + Block.Shelf + Block.RowSplit"

emm_AGB_by_geno <- emmeans(specs = c("GN"), type = "response",
                                    FE_analysis(my_data %>% mutate("GN" = interaction(Genotype, Nema)),
                                                formula = "Above.ground.biomass ~ GN",  plot = FALSE, block = FALSE,
                                                singular.ok = TRUE, family = lognormal,
                                                return_model = TRUE)) %>%
  as.data.frame() %>% rowwise() %>%
  mutate(Genotype = unlist(strsplit(split = "\\.", as.character(GN)))[1],
         Nema = unlist(strsplit(split = "\\.", as.character(GN)))[2])

p_AGB_geno <- ggplot(emm_AGB_by_geno) + 
  aes(x = Nema, y = response, ymax = response + SE, ymin = response - SE) +
  aes(fill = Genotype) + geom_pointrange(position = position_dodge(width = 0.33)) + 
  # geom_text(data = . %>% filter(Nema == "+"), hjust = -1,
  #           x = "+", nudge_x = 0.5, aes(label = Genotype)) + 
  # geom_text(data = . %>% filter(Nema == "-"), hjust = 2,
  #           x = "-", nudge_x = -0.5, aes(label = Genotype)) + 
  geom_line(aes(group = Genotype), position = position_dodge(width = 0.33)) +
  scale_x_discrete(labels = c("-" = "Uninfected", "+" = "Infected")) +
  xlab("Nematode infection status") +
  ylab("Above.ground.biomass") +
  theme_minimal() + ylim(c(0,0.25)) +
  guides(fill = guide_legend(ncol = 2)) +
  ggtitle("A") + theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot")


emm_AGB_by_rhizo <- emmeans(specs = c("RN"), type = "response",
                           FE_analysis(my_data %>% mutate("RN" = interaction(Rhizo, Nema)),
                                       formula = "Above.ground.biomass ~ RN", plot = FALSE, block = FALSE,
                                       singular.ok = TRUE, family = lognormal,
                                       return_model = TRUE)) %>%
  as.data.frame() %>% rowwise() %>%
  mutate(Rhizo = unlist(strsplit(split = "\\.", as.character(RN)))[1],
         Nema = unlist(strsplit(split = "\\.", as.character(RN)))[2])

p_AGB_rhizo <- ggplot(emm_AGB_by_rhizo) + 
  aes(x = Nema, y = response, ymax = response + SE, ymin = response - SE) +
  aes(fill = Rhizo) + geom_pointrange(shape = 21, position = position_dodge(width = 0.33)) + 
  geom_line(aes(group = Rhizo), position = position_dodge(width = 0.33)) +
  scale_x_discrete(labels = c("-" = "Uninfected", "+" = "Infected")) +
  xlab("Nematode infection status") +
  ylab("Above ground biomass") +
  theme_minimal() + ylim(c(0,0.25)) +
  guides(fill = guide_legend(ncol = 2)) +
  ggtitle("B") + ylim(c(0,0.25)) +
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot")


p_AGB_variance <- ggplot(AGB_RE %>% filter(`Chi Df` == 1)%>%
                           filter(Effect != "(Intercept)" & Effect != "Block.Position1" & Effect != "scale(Volume.mm3)")) +
  aes(fill = Effect, y = Variance, x = `Chi Df`, label = sig) +
  geom_bar(position = "fill", stat="identity") +
  geom_text(data = . %>%
              filter(!is.na(Variance), Effect != "Rhizo") %>% 
              mutate(pct = Variance / sum(Variance, na.rm = TRUE)),
            aes(label = paste(percent(pct), ifelse(sig != ".", sig, ""))),
            position = position_fill(vjust =0.5)) +
  theme_minimal() +
  ggtitle("C") + theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())

emt_AGB_sev <- emmip(FE_analysis(my_data %>% filter(Nema == "+"), plot = FALSE,
                                 formula = "Above.ground.biomass ~ Galls*Rhizo + Genotype*Galls", block = FALSE,
                                 family_distribution = lognormal, singular.ok = TRUE, return_model = TRUE),
  ~Galls:Rhizo, var = "Galls", specs = c("Rhizo", "Galls"), type = "response", at = list(Galls = c(seq(1,61,2.5))), CIs = TRUE, plotit = FALSE, rg.limit = 1000000) %>% as.data.frame()
  
p_AGB_sev <- ggplot(emt_AGB_sev) + 
  geom_line(aes(x=Galls, y=yvar, color=Rhizo)) +
  ylab("Above ground biomass (g)") +
  geom_ribbon(aes(x = Galls, ymax=UCL, ymin=LCL, fill=Rhizo), alpha = 0.2) +
  theme_minimal() +
  ggtitle("D") +
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot")

p_AGB <- grid.arrange(p_AGB_geno, p_AGB_rhizo, p_AGB_variance, p_AGB_sev,
                      ncol = 3, nrow=2, layout_matrix = rbind(c(1,1,2),c(3,4,4)), widths = c(2,1,3))
