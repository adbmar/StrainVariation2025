################################################################################
# Author: Addison Dale Buxton-Martin
# 
# Description: This script runs the analyses associated with the publication ___
# This analysis will rely on other scripts in its encasing directory to import
# data and define helper functions and then perform analyses.
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
if(!exists("my_data")){source(file.path(dir_main,"LoadData.R"))}
if(!exists("FE_analysis")){source(file.path(dir_main,"AnalysisFunctions.R"))}

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


block_formula <- "+ Block.Position1 + Block.Shelf + Block.RowSplit"

###########################################################
### EFFECT OF MUTUALIST GENETIC VARIATION ON PARASITISM ###
###########################################################

#########
### Does rhizobia variation impact host parasite load?
#########

Gall_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                             formula = "Galls ~ Genotype + Rhizo",
                             family_distribution = nbinom2, singular.ok = TRUE)

Gall_vc_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                          formula = "Galls ~ Genotype + Rhizo + scale(Volume.mm3)",
                          family_distribution = nbinom2, singular.ok = TRUE)
Gall_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"),
                       family_distribution = nbinom1,
                       fullest_formula = "Galls ~ (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo)",
                       formula_table =
                         rbind(
                           c("Effect" = "Genotype:Rhizo",
                             "Full" = "Galls ~ (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo)",
                             "Reduced" = "Galls ~ (1|Genotype) + (1|Rhizo)"),
                           c("Effect" = "Rhizo", 
                             "Full" = "Galls ~ (1|Genotype) + (1|Rhizo)",
                             "Reduced" = "Galls ~ (1|Genotype)"),
                           c("Effect" = "Genotype",
                             "Full" = "Galls ~ (1|Genotype) + (1|Rhizo)",
                             "Reduced" = "Galls ~ (1|Rhizo)")))

Gall_vc_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"),
                          formula_table = 
                            rbind(
                              c("Effect" = "Genotype:Rhizo",
                                "Full" = "Galls ~ scale(Volume.mm3) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo)",
                                "Reduced" = "Galls ~ scale(Volume.mm3) + (1|Genotype) + (1|Rhizo)"),
                              c("Effect" = "Rhizo", 
                                "Full" = "Galls ~ scale(Volume.mm3) + (1|Genotype) + (1|Rhizo)",
                                "Reduced" = "Galls ~ scale(Volume.mm3) + (1|Genotype)"),
                              c("Effect" = "Genotype",
                                "Full" = "Galls ~ scale(Volume.mm3) + (1|Genotype) + (1|Rhizo)",
                                "Reduced" = "Galls ~ scale(Volume.mm3) + (1|Rhizo)")),
                          fullest_formula = "Galls ~ scale(Volume.mm3) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo)",
                          family_distribution = nbinom1)

Gall_sub1_FE <- FE_analysis(df = subcluster_1, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo")
Gall_sub2_FE <- FE_analysis(df = subcluster_2, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo")
Gall_sub3_FE <- FE_analysis(df = subcluster_3, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo")
Gall_sub4_FE <- FE_analysis(df = subcluster_4, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo")
Gall_sub5_FE <- FE_analysis(df = subcluster_5, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo")
Gall_sub1_vc_FE <- FE_analysis(df = subcluster_1, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo + scale(Volume.mm3)")
Gall_sub2_vc_FE <- FE_analysis(df = subcluster_2, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo + scale(Volume.mm3)")
Gall_sub3_vc_FE <- FE_analysis(df = subcluster_3, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo + scale(Volume.mm3)")
Gall_sub4_vc_FE <- FE_analysis(df = subcluster_4, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo + scale(Volume.mm3)")
Gall_sub5_vc_FE <- FE_analysis(df = subcluster_5, family_distribution = nbinom1, singular.ok = TRUE,
                            formula = "Galls ~ Genotype * Rhizo + scale(Volume.mm3)")

###########################################################
### EFFECT OF MUTUALIST GENETIC VARIATION ON PARASITISM ###
###########################################################

#########
### Does rhizobial variation impact the mutualism's response to parasite infection?
#########

Nods_FE <- FE_analysis(df = my_data,
                       formula = "Total.Nodules ~ Genotype*Nema + Rhizo*Nema",
                       family_distribution = nbinom1, singular.ok = TRUE)
Nods_FE_noblock <- FE_analysis(df = my_data, block = FALSE,
                       formula = "Total.Nodules ~ Genotype*Nema + Rhizo*Nema",
                       family_distribution = nbinom1, singular.ok = TRUE)

Nods_vc_FE <- FE_analysis(df = my_data,
                       formula = "Total.Nodules ~ Genotype*Nema + Rhizo*Nema + scale(Volume.mm3)",
                       family_distribution = nbinom1, singular.ok = TRUE)
Nods_vc_FE_noblock <- FE_analysis(df = my_data, block = FALSE,
                          formula = "Total.Nodules ~ Genotype*Nema + Rhizo*Nema + scale(Volume.mm3)",
                          family_distribution = nbinom1, singular.ok = TRUE)
Nods_RE <- RE_analysis(df = my_data,
                       fullest_formula = "Total.Nodules ~ Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo) + (1|Nema:Genotype:Rhizo)",
                       formula_table = generate_formula_table("Total.Nodules"))
Nods_vc_RE <- RE_analysis(df = my_data, family_distribution = nbinom1, 
                            fullest_formula = "Total.Nodules ~ scale(Volume.mm3) + Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo) + (1|Nema:Genotype:Rhizo)",
                            formula_table = generate_formula_table("Total.Nodules", volume = "scale"))


Nods_sub1_FE <- FE_analysis(df = subcluster_1, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Nema")
Nods_sub2_FE <- FE_analysis(df = subcluster_2, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Nema")
Nods_sub3_FE <- FE_analysis(df = subcluster_3, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Nema")
Nods_sub4_FE <- FE_analysis(df = subcluster_4, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Nema")
Nods_sub5_FE <- FE_analysis(df = subcluster_5, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Nema")
Nods_sub1_vc_FE <- FE_analysis(df = subcluster_1, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Nema + scale(Volume.mm3)")
Nods_sub2_vc_FE <- FE_analysis(df = subcluster_2, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                               formula = "Total.Nodules ~ Genotype * Rhizo * Nema + scale(Volume.mm3)")
Nods_sub3_vc_FE <- FE_analysis(df = subcluster_3, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                               formula = "Total.Nodules ~ Genotype * Rhizo * Nema + scale(Volume.mm3)")
Nods_sub4_vc_FE <- FE_analysis(df = subcluster_4, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                               formula = "Total.Nodules ~ Genotype * Rhizo * Nema + scale(Volume.mm3)")
Nods_sub5_vc_FE <- FE_analysis(df = subcluster_5, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                               formula = "Total.Nodules ~ Genotype * Rhizo * Nema + scale(Volume.mm3)")


#########
### Does rhizobia variation impact the mutualism's response to host parasite load?
#########
Nods_sev_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                       formula = "Total.Nodules ~ Genotype*scale(Galls) + Rhizo*scale(Galls)",
                       family_distribution = nbinom1, singular.ok = TRUE)

Nods_sev_vc_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                          formula = "Total.Nodules ~ Genotype*scale(Galls) + Rhizo*scale(Galls) + scale(Volume.mm3)",
                          family_distribution = nbinom1, singular.ok = TRUE)

Nods_sev_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = gaussian,
                           fullest_formula = "scale(Total.Nodules) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)",
                           formula_table = generate_formula_table("Total.Nodules", galls = "scale", mod_trait = "scale"))

Nods_sev_vc_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = gaussian,
                                fullest_formula = "scale(Total.Nodules) ~ scale(Volume.mm3) + scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)",
                                formula_table = generate_formula_table("Total.Nodules", volume = "scale", galls = "scale", mod_trait = "scale"))


Nods_sev_sub1_FE <- FE_analysis(df = subcluster_1, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Galls")
Nods_sev_sub2_FE <- FE_analysis(df = subcluster_2, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Galls")
Nods_sev_sub3_FE <- FE_analysis(df = subcluster_3, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Galls")
Nods_sev_sub4_FE <- FE_analysis(df = subcluster_4, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Galls")
Nods_sev_sub5_FE <- FE_analysis(df = subcluster_5, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                            formula = "Total.Nodules ~ Genotype * Rhizo * Galls")
Nods_sev_sub1_vc_FE <- FE_analysis(df = subcluster_1, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                               formula = "Total.Nodules ~ Genotype * Rhizo * Galls + scale(Volume.mm3)")
Nods_sev_sub2_vc_FE <- FE_analysis(df = subcluster_2, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                                   formula = "Total.Nodules ~ Genotype * Rhizo * Galls + scale(Volume.mm3)")
Nods_sev_sub3_vc_FE <- FE_analysis(df = subcluster_3, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                                   formula = "Total.Nodules ~ Genotype * Rhizo * Galls + scale(Volume.mm3)")
Nods_sev_sub4_vc_FE <- FE_analysis(df = subcluster_4, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                                   formula = "Total.Nodules ~ Genotype * Rhizo * Galls + scale(Volume.mm3)")
Nods_sev_sub5_vc_FE <- FE_analysis(df = subcluster_5, family_distribution = nbinom1, singular.ok = TRUE, block = FALSE,
                                   formula = "Total.Nodules ~ Genotype * Rhizo * Galls + scale(Volume.mm3)")


#########
### Does rhizobia variation impact the host's fitness conseqeunces for parasite infection?
#########
AGB_FE <- FE_analysis(df = my_data,
                           formula = "Above.ground.biomass ~ Genotype*Nema + Rhizo*Nema",
                           family_distribution = lognormal, singular.ok = TRUE)

AGB_RE <- RE_analysis(df = my_data, family_distribution = lognormal,
                      fullest_formula = "Above.ground.biomass ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Genotype:Rhizo)",
                      formula_table = generate_formula_table("Above.ground.biomass"))

AGB_nods_FE <- FE_analysis(df = my_data,
                           formula = "Above.ground.biomass ~ Nema*Rhizo*Total.Nodules",
                           family_distribution = lognormal, singular.ok = TRUE)
AGB_nods_RE <- RE_analysis(df = my_data, family_distribution = lognormal, singular.ok = TRUE,
                           fullest_formula = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (1|Nema:Rhizo) + (0+Total.Nodules|Rhizo) + (0+Total.Nodules|Nema) + (0+Total.Nodules|Nema:Rhizo)",
                           formula_table =
                             rbind(
                               c("Effect" = "Nema.Rhizo.1",
                                 "Full" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (1|Nema:Rhizo) + (0+Total.Nodules|Rhizo) + (0+Total.Nodules|Nema) + (0+Total.Nodules|Nema:Rhizo)",
                                 "Reduced" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (1|Nema:Rhizo) + (0+Total.Nodules|Rhizo) + (0+Total.Nodules|Nema)"),
                               c("Effect" = "Rhizo.1",
                                 "Full" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (1|Nema:Rhizo) + (0+Total.Nodules|Rhizo) + (0+Total.Nodules|Nema)",
                                 "Reduced" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (1|Nema:Rhizo) + (0+Total.Nodules|Nema)"),
                                 c("Effect" = "Nema.1",
                                   "Full" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (1|Nema:Rhizo) + (0+Total.Nodules|Rhizo) + (0+Total.Nodules|Nema)",
                                   "Reduced" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (1|Nema:Rhizo) + (0+Total.Nodules|Rhizo)"),
                               c("Effect" = "Nema.Rhizo",
                                 "Full" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (1|Nema:Rhizo) + (0+Total.Nodules|Rhizo) + (0+Total.Nodules|Nema)",
                                 "Reduced" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo) + (0+Total.Nodules|Rhizo) + (0+Total.Nodules|Nema)"),
                               c("Effect" = "Rhizo",
                                 "Full" = "Above.ground.biomass ~ Nema + Total.Nodules + (1|Rhizo)",
                                 "Reduced" = "Above.ground.biomass ~ Nema + Total.Nodules")))


AGB_sev_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                      formula = "Above.ground.biomass ~ scale(Galls)*Genotype + scale(Galls)*Rhizo",
                      family_distribution = lognormal, singular.ok = TRUE)

AGB_sev_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = lognormal,
                      fullest_formula = "Above.ground.biomass ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)",
                      formula_table = generate_formula_table("Above.ground.biomass", galls = "scale"))


### Double checking the pvalue of 1 in the above model results
full_test_formula <- Above.ground.biomass ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + Block.Position1 + Block.Shelf + Block.RowSplit
reduced_test_formula <- Above.ground.biomass ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + Block.Position1 + Block.Shelf + Block.RowSplit
simulateResiduals(plot = T, 
                  glmmTMB(my_data %>% filter(Nema == "+"), family = gaussian,
                          formula = full_test_formula))
anova(glmmTMB(my_data %>% filter(Nema == "+"), family = gaussian, formula = full_test_formula),
      glmmTMB(my_data %>% filter(Nema == "+"), family = gaussian, formula = reduced_test_formula), 
      test = "LRT")

glmmTMB(my_data %>% filter(Nema == "+"), family = lognormal, formula = reduced_test_formula)
glmmTMB(my_data %>% filter(Nema == "+"), family = lognormal, formula = full_test_formula)


SRrat_FE <- FE_analysis(df = my_data,
                        formula = "SvR_ratio ~ Genotype*Nema + Rhizo*Nema",
                        family_distribution = lognormal, singular.ok = TRUE)

SRrat_FE_noblock <- FE_analysis(df = my_data, block = FALSE,
                        formula = "SvR_ratio ~ Genotype*Nema + Rhizo*Nema",
                        family_distribution = lognormal, singular.ok = TRUE)

SRrat_sev_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                        formula = "SvR_ratio ~ scale(Galls)*Genotype + scale(Galls)*Rhizo",
                        family_distribution = lognormal, singular.ok = TRUE)

SRrat_RE <- RE_analysis(df = my_data %>% filter(!is.na(SvR_ratio) & SvR_ratio>0), family_distribution = lognormal,
                          fullest_formula = "SvR_ratio ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Genotype:Rhizo)",
                          formula_table = generate_formula_table("SvR_ratio"))

### Testing a subtest that had model convergence problems in the above
full_test_formula <- SvR_ratio ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)
reduced_test_formula <- SvR_ratio ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Nema)
simulateResiduals(plot = T, 
                  glmmTMB(my_data, family = gaussian(link = "log"),
                          formula = full_test_formula))
anova(glmmTMB(my_data, family = gaussian(link = "log"), formula = full_test_formula),
      glmmTMB(my_data, family = gaussian(link = "log"), formula = reduced_test_formula), 
      test = "LRT")
#results show a very high p value even here, so safe to assume the pvalue is 1 in the RE_analysis test


SRrat_sev_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = gaussian(link = "identity"),
                              fullest_formula = "scale(SvR_ratio) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)",
                              formula_table = generate_formula_table("SvR_ratio", galls = "scale", mod_trait = "scale"))

roottip_FE <- FE_analysis(df = my_data,
                           formula = "Number.of.Root.Tips ~ Nema*Genotype + Nema*Rhizo",
                           family_distribution = nbinom1, singular.ok = TRUE)

exploratory <- FE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = nbinom2,
                          formula = "Galls ~ Genotype*Number.of.Root.Tips + Rhizo*Number.of.Root.Tips",
                          singular.ok = TRUE)
ggplot(
  my_data %>% group_by(Nema, Genotype, Rhizo) %>% filter(Nema == "+") %>%
    summarize(galls_mean = mean(Galls), galls_sd = sd(Galls), 
              root_mean = mean(Number.of.Root.Tips), root_sd = sd(Number.of.Root.Tips))) +
  aes(x = galls_mean, xmin = galls_mean - galls_sd, xmax = galls_mean + galls_sd,
      y = root_mean, ymin = root_mean - root_sd, ymax = root_mean + root_sd,
      fill = Rhizo, color = Rhizo, shape = Nema) +
  geom_errorbarh() + geom_pointrange() +
  scale_shape_manual(values = c("+" = 21, "-" = 25)) +
  theme_classic()

ggplot(
  my_data %>% group_by(Nema, Rhizo) %>% filter(Nema == "+") %>%
    summarize(galls_mean = mean(Galls), galls_sd = sd(Galls), 
              root_mean = mean(Number.of.Root.Tips), root_sd = sd(Number.of.Root.Tips))) +
  aes(x = galls_mean, xmin = galls_mean - galls_sd, xmax = galls_mean + galls_sd,
      y = root_mean, ymin = root_mean - root_sd, ymax = root_mean + root_sd,
      fill = Rhizo, shape = Nema) +
  geom_errorbarh() + geom_pointrange() +
  scale_shape_manual(values = c("+" = 21, "-" = 25)) +
  theme_classic()

roottip_FE_noblock <- FE_analysis(df = my_data, block=FALSE,
                          formula = "Number.of.Root.Tips ~ Nema*Genotype + Nema*Rhizo",
                          family_distribution = nbinom1, singular.ok = TRUE)

roottip_sev_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                          formula = "Number.of.Root.Tips ~ scale(Galls)*Genotype + scale(Galls)*Rhizo",
                          family_distribution = nbinom1, singular.ok = TRUE)

roottip_RE <- RE_analysis(df = my_data, family_distribution = nbinom1,
                            fullest_formula = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Genotype:Rhizo)",
                            formula_table = generate_formula_table("Number.of.Root.Tips"))

full_test_formula <- Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype)
reduced_test_formula <- Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo)
simulateResiduals(plot = T, 
                  glmmTMB(my_data, family = nbinom1,
                          formula = full_test_formula))
anova(glmmTMB(my_data, family = nbinom1, formula = full_test_formula),
      glmmTMB(my_data, family = nbinom1, formula = reduced_test_formula), 
      test = "LRT")

roottip_RE_2 <- RE_analysis(df = my_data, family_distribution = nbinom1,
                          fullest_formula = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Genotype:Rhizo)",
                          formula_table = 
                              rbind(
                                c("Effect" = "Nema:Genotype:Rhizo",
                                  "Full" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Genotype:Rhizo)",
                                  "Reduced" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype)"),
                                c("Effect" = "Nema:Genotype",
                                  "Full" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype)",
                                  "Reduced" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo)"),
                                c("Effect" = "Nema:Rhizo",
                                  "Full" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype)",
                                  "Reduced" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Genotype:Rhizo) + (1|Nema:Genotype)"),
                                c("Effect" = "Genotype:Rhizo",
                                  "Full" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype)",
                                  "Reduced" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Nema:Genotype)"),
                                c("Effect" = "Genotype",
                                  "Full" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype)",
                                  "Reduced" = "Number.of.Root.Tips ~ Nema + (1|Rhizo)"),
                                c("Effect" = "Rhizo",
                                  "Full" = "Number.of.Root.Tips ~ Nema + (1|Rhizo) + (1|Genotype)",
                                  "Reduced" = "Number.of.Root.Tips ~ Nema + (1|Genotype)")))
roottip_sev_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = gaussian(link = "identity"),
                              fullest_formula = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)",
                              formula_table = 
                                rbind(
                                  c("Effect" = "Genotype.Rhizo.1",
                                    "Full" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)",
                                    "Reduced" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)"),
                                  c("Effect" = "Genotype.1",
                                    "Full" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                    "Reduced" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Rhizo)"),
                                  c("Effect" = "Rhizo.1",
                                    "Full" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                    "Reduced" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype)"),
                                  c("Effect" = "Genotype.Rhizo",
                                    "Full" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                    "Reduced" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)"),
                                  c("Effect" = "Genotype",
                                    "Full" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo)",
                                    "Reduced" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Rhizo)"),
                                  c("Effect" = "Rhizo",
                                    "Full" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype) + (1|Rhizo)",
                                    "Reduced" = "scale(Number.of.Root.Tips) ~ scale(Galls) + (1|Genotype)")))

roottip_sev_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = nbinom2,
                              fullest_formula = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo)",
                              formula_table = 
                                rbind(
                                  c("Effect" = "Genotype.Rhizo.1",
                                    "Full" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo)",
                                    "Reduced" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)"),
                                  c("Effect" = "Genotype.1",
                                    "Full" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
                                    "Reduced" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Rhizo)"),
                                  c("Effect" = "Rhizo.1",
                                    "Full" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
                                    "Reduced" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype)"),
                                  c("Effect" = "Genotype.Rhizo",
                                    "Full" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
                                    "Reduced" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)"),
                                  c("Effect" = "Genotype",
                                    "Full" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo)",
                                    "Reduced" = "Number.of.Root.Tips ~ Galls + (1|Rhizo)"),
                                  c("Effect" = "Rhizo",
                                    "Full" = "Number.of.Root.Tips ~ Galls + (1|Genotype) + (1|Rhizo)",
                                    "Reduced" = "Number.of.Root.Tips ~ Galls + (1|Genotype)")))
### Testing a subtest that had model convergence problems in the above
full_test_formula <- Number.of.Branch.Points ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)
reduced_test_formula <- Number.of.Branch.Points ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)
simulateResiduals(plot = T, 
                  glmmTMB(my_data %>% filter(Nema == "+"), family = gaussian(link = "log"),
                          formula = full_test_formula))
anova(glmmTMB(my_data, family = nbinom1, formula = full_test_formula),
      glmmTMB(my_data, family = nbinom1, formula = reduced_test_formula), 
      test = "LRT")
#results show a very high p value even here, so safe to assume the pvalue is 1 in the RE_analysis test

branchfreq_FE <- FE_analysis(df = my_data,
                              formula = "scale(Branching.frequency.per.mm) ~ Nema*Genotype + Nema*Rhizo",
                             family_distribution = gaussian, singular.ok = TRUE)
branchfreq_FE_noblock <- FE_analysis(df = my_data, block = FALSE,
                             formula = "scale(Branching.frequency.per.mm) ~ Nema*Genotype + Nema*Rhizo",
                             family_distribution = gaussian, singular.ok = TRUE)
branchfreq_RE <- RE_analysis(df = my_data, family_distribution = gaussian,
                             fullest_formula = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema) + (1|Nema:Genotype:Rhizo)",
                             formula_table = 
                               rbind(
                                 c("Effect" = "Nema:Genotype:Rhizo",
                                   "Full" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema) + (1|Nema:Genotype:Rhizo)",
                                   "Reduced" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)"),
                                 c("Effect" = "Genotype:Nema",
                                   "Full" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)",
                                   "Reduced" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo)"),
                                 c("Effect" = "Nema:Rhizo",
                                   "Full" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)",
                                   "Reduced" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Genotype:Rhizo) + (1|Genotype:Nema)"),
                                 c("Effect" = "Genotype:Rhizo",
                                   "Full" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)",
                                   "Reduced" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Nema)"),
                                 c("Effect" = "Genotype",
                                   "Full" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype)",
                                   "Reduced" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo)"),
                                 c("Effect" = "Rhizo",
                                   "Full" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype)",
                                   "Reduced" = "scale(Branching.frequency.per.mm) ~ Nema + (1|Genotype)")))
full_test_formula <- scale(Branching.frequency.per.mm) ~ Nema + (1|Rhizo) + (1|Genotype) + Block.Position1 + Block.Shelf + Block.RowSplit
reduced_test_formula <- scale(Branching.frequency.per.mm) ~ Nema + (1|Genotype) + Block.Position1 + Block.Shelf + Block.RowSplit
simulateResiduals(plot = T, 
                  glmmTMB(my_data, family = gaussian,
                          formula = full_test_formula))
anova(glmmTMB(my_data, family = gaussian, formula = full_test_formula),
      glmmTMB(my_data, family = gaussian, formula = reduced_test_formula), 
      test = "LRT")

branchfreq_sev_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                             formula = "scale(Branching.frequency.per.mm) ~ scale(Galls)*Genotype + scale(Galls)*Rhizo",
                             family_distribution = gaussian, singular.ok = TRUE)
branchfreq_sev_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = gaussian,
                              fullest_formula = "scale(Branching.frequency.per.mm) ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo)",
                              formula_table = 
                                rbind(
                                  c("Effect" = "Genotype.Rhizo.1",
                                    "Full" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)",
                                    "Reduced" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)"),
                                  c("Effect" = "Genotype.1",
                                    "Full" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                    "Reduced" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Rhizo)"),
                                  c("Effect" = "Rhizo.1",
                                    "Full" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                    "Reduced" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype)"),
                                  c("Effect" = "Genotype.Rhizo",
                                    "Full" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                    "Reduced" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)"),
                                  c("Effect" = "Genotype",
                                    "Full" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo)",
                                    "Reduced" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Rhizo)"),
                                  c("Effect" = "Rhizo",
                                    "Full" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype) + (1|Rhizo)",
                                    "Reduced" = "scale(Branching.frequency.per.mm) ~ scale(Galls) + (1|Genotype)")))

branchnum_FE <- FE_analysis(df = my_data,
                             formula = "Number.of.Branch.Points ~ Nema*Genotype + Nema*Rhizo",
                             family_distribution = nbinom1, singular.ok = TRUE)
branchnum_FE_noblock <- FE_analysis(df = my_data, block = FALSE,
                            formula = "Number.of.Branch.Points ~ Nema*Genotype + Nema*Rhizo",
                            family_distribution = nbinom1, singular.ok = TRUE)

branchnum_RE <- RE_analysis(df = my_data, family_distribution = nbinom1,
                            fullest_formula = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema) + (1|Nema:Genotype:Rhizo)",
                            formula_table = 
                              rbind(
                                c("Effect" = "Nema:Genotype:Rhizo",
                                  "Full" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema) + (1|Nema:Genotype:Rhizo)",
                                  "Reduced" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)"),
                                c("Effect" = "Genotype:Nema",
                                  "Full" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)",
                                  "Reduced" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo)"),
                                c("Effect" = "Nema:Rhizo",
                                  "Full" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)",
                                  "Reduced" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Genotype:Rhizo) + (1|Genotype:Nema)"),
                                c("Effect" = "Genotype:Rhizo",
                                  "Full" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Rhizo) + (1|Genotype:Nema)",
                                  "Reduced" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype) + (1|Nema:Rhizo) + (1|Genotype:Nema)"),
                                c("Effect" = "Genotype",
                                  "Full" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype)",
                                  "Reduced" = "Number.of.Branch.Points ~ Nema + (1|Rhizo)"),
                                c("Effect" = "Rhizo",
                                  "Full" = "Number.of.Branch.Points ~ Nema + (1|Rhizo) + (1|Genotype)",
                                  "Reduced" = "Number.of.Branch.Points ~ Nema + (1|Genotype)")))

branchnum_sev_FE <- FE_analysis(df = my_data %>% filter(Nema == "+"),
                            formula = "Number.of.Branch.Points ~ scale(Galls)*Genotype + scale(Galls)*Rhizo",
                            family_distribution = nbinom1, singular.ok = TRUE)
branchnum_sev_RE <- RE_analysis(df = my_data %>% filter(Nema == "+"), family_distribution = gaussian,
                                 fullest_formula = "scale(Number.of.Branch.Points) ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo)",
                                 formula_table = 
                                   rbind(
                                     c("Effect" = "Genotype.Rhizo.1",
                                       "Full" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo) + (0+scale(Galls)|Genotype:Rhizo)",
                                       "Reduced" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)"),
                                     c("Effect" = "Genotype.1",
                                       "Full" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                       "Reduced" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Rhizo)"),
                                     c("Effect" = "Rhizo.1",
                                       "Full" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                       "Reduced" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype)"),
                                     c("Effect" = "Genotype.Rhizo",
                                       "Full" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)",
                                       "Reduced" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo) + (0+scale(Galls)|Genotype) + (0+scale(Galls)|Rhizo)"),
                                     c("Effect" = "Genotype",
                                       "Full" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo)",
                                       "Reduced" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Rhizo)"),
                                     c("Effect" = "Rhizo",
                                       "Full" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype) + (1|Rhizo)",
                                       "Reduced" = "scale(Number.of.Branch.Points) ~ scale(Galls) + (1|Genotype)")))
