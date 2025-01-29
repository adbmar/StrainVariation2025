################################################################################
# Author: Addison Dale Buxton-Martin
# 
# Description: This script is part of a broader group of scripts
# that documents the analysis associated with the publication found at:

# This script loads data used for an analysis of data for a
# publication. This will import data from a csv included in its encasing
# directory and perform preliminary data checks and analyses. This script
# is also called by other scripts that require the data being already loaded.
#
################################################################################

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
  "tidyverse","corrplot"))


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


######################
### Importing data ###
######################
data <- read.csv(file.path(dir_main, "cleaned_data_mod.csv"))
data <- data %>% mutate(Block.Position1 = substring(Block.Position, 1, 1), #Changing Block.Position to avoid aliased coefficients
                        Block.Position2 = substring(Block.Position, 2, 2), #Changing Block.Position to avoid aliased coefficients
                        Block.Split = if_else(is.na(Block.Col) | is.na(Block.Row), NA,
                                              case_when(Block.Col %in% c(2,3) & Block.Row %in% c(2,3) ~ "A",
                                                        Block.Col %in% c(5,6,7) & Block.Row %in% c(4,5) ~ "B",
                                                        Block.Col %in% c(8,9,10) & Block.Row %in% c(6,7) ~ "C",
                                                        Block.Col %in% c(12,13) & Block.Row %in% c(2,3) ~ "D",
                                                        Block.Col %in% c(2,3) & Block.Row %in% c(4,5) ~ "E",
                                                        Block.Col %in% c(5,6,7) & Block.Row %in% c(6,7) ~ "F",
                                                        Block.Col %in% c(8,9,10) & Block.Row %in% c(2,3) ~ "G",
                                                        Block.Col %in% c(12,13) & Block.Row %in% c(4,5) ~ "H",
                                                        Block.Col %in% c(2,3) & Block.Row %in% c(6,7) ~ "I",
                                                        Block.Col %in% c(5,6,7) & Block.Row %in% c(2,3) ~ "J",
                                                        Block.Col %in% c(8,9,10) & Block.Row %in% c(4,5) ~ "K",
                                                        Block.Col %in% c(12,13) & Block.Row %in% c(6,7) ~ "L")))
data <- data %>% 
  mutate(Block.Large = interaction(Block.Shelf, Block.Position1)) %>%
  mutate(Block.Small = Block.Split) %>%
  mutate(Total.Nodules = as.integer(Nodules) + as.integer(Nods.Removed)) %>%
  mutate(Above.ground.biomass = as.numeric(Above.ground.biomass)) %>%
  mutate(Block.ColSplit = case_when(Block.Col %in% c(2,3) ~ "A",
                                    Block.Col %in% c(5,6,7) ~ "B",
                                    Block.Col %in% c(8,9,10) ~ "C",
                                    Block.Col %in% c(12,13) ~ "D",
                                    is.na(Block.Col) ~ NA,
                                    .default = "E")) %>%
  mutate(Block.RowSplit = case_when(Block.Row %in% c(2,3) ~ "A",
                                    Block.Row %in% c(4,5) ~ "B",
                                    Block.Row %in% c(6,7) ~ "C",
                                    .default = "D")) %>%
  mutate(Block.Block = interaction(Block.Position1, Block.ColSplit)) %>%
  mutate(Block.Block2 = interaction(Block.Shelf, Block.Position1, Block.ColSplit))

RV_data_root_nema_neg <- read.csv(file.path(dir_main, "RhizoVision_RootNemaNeg.csv"))
RV_data_root_nema_pos <- read.csv(file.path(dir_main, "RhizoVision_RootNemaPos.csv"))

RV_data_root_nema_neg <- RV_data_root_nema_neg %>% 
  select(File.Name, Number.of.Root.Tips, Number.of.Branch.Points, Total.Root.Length.mm,
         Branching.frequency.per.mm, Network.Area.mm2, Average.Diameter.mm,
         Median.Diameter.mm, Maximum.Diameter.mm, Perimeter.mm, Volume.mm3,
         Surface.Area.mm2) %>%
  rename("ids" = File.Name)

RV_data_root_nema_pos <- RV_data_root_nema_pos %>%
  select(File.Name, Number.of.Root.Tips, Number.of.Branch.Points, Total.Root.Length.mm,
         Branching.frequency.per.mm, Network.Area.mm2, Average.Diameter.mm,
         Median.Diameter.mm, Maximum.Diameter.mm, Perimeter.mm, Volume.mm3,
         Surface.Area.mm2) %>%
  rename("ids" = File.Name)

RV_data <- rbind(RV_data_root_nema_pos, RV_data_root_nema_neg, by = "ids") %>% filter(ids != "ids")
RV_data$ids <- gsub(".png", "", RV_data$ids)

data <- full_join(RV_data, data, by = "ids")



data$Above.ground.biomass <- as.numeric(data$Above.ground.biomass)
data$Number.of.Root.Tips <- as.numeric(data$Number.of.Root.Tips)
data$Number.of.Branch.Points <- as.numeric(data$Number.of.Branch.Points)
data$Total.Root.Length.mm <- as.numeric(data$Total.Root.Length.mm)
data$Branching.frequency.per.mm <- as.numeric(data$Branching.frequency.per.mm)
data$Network.Area.mm2 <- as.numeric(data$Network.Area.mm2)
data$Network.Area.mm2 <- as.numeric(data$Network.Area.mm2)
data$Average.Diameter.mm <- as.numeric(data$Average.Diameter.mm)
data$Median.Diameter.mm <- as.numeric(data$Median.Diameter.mm)
data$Maximum.Diameter.mm <- as.numeric(data$Maximum.Diameter.mm)
data$Perimeter.mm <- as.numeric(data$Perimeter.mm)
data$Volume.mm3 <- as.numeric(data$Volume.mm3)
data$Surface.Area.mm2 <- as.numeric(data$Surface.Area.mm2)
data$Nema <- as.factor(data$Nema)
data$Rhizo <- as.factor(data$Rhizo)
data$Genotype <- as.factor(data$Genotype)
data$Block.Small <- as.factor(data$Block.Small)
data$Block.Large <- as.factor(data$Block.Large)
data$Block.Shelf <- as.factor(data$Block.Shelf)
data$Block.Side <- as.factor(data$Block.Side)
data$Block.Position <- as.factor(data$Block.Position)
data$Block.Col <- as.factor(data$Block.Col)
data$Block.ColSplit <- as.factor(data$Block.ColSplit)
data$Block.Block <- as.factor(data$Block.Block)

data <- data %>% mutate(SvR_ratio = as.numeric(Above.ground.biomass/Volume.mm3))

### Single data point filtering
#sanity check of volume and surface area
ggplot(data) + aes(x = Surface.Area.mm2, y = Volume.mm3, label = ids) + geom_point() + theme_classic()
ggplot(data) + aes(x = Surface.Area.mm2, y = Volume.mm3, label = ids) + geom_point() + theme_classic() + geom_text()

#Overall checking that data looks okay
ggplot(
  data %>% mutate(across(where(is.numeric), scale)) %>% pivot_longer(cols = data %>% select(where(is.numeric)) %>% colnames())
) + aes(x = name, y = value, label = ids) + geom_point() + 
  geom_text(data = . %>% filter(value >= 3 | value <= -3), position = position_nudge(x = 0.6)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90))

#filter out weird surface area volume points
data <- data %>% filter(!ids %in% c("id490", "id487", "id065", "id041", "id065", "id041", "id060", "id562"))

ggplot(data) + aes(x = Surface.Area.mm2, y = Volume.mm3, label = ids) + geom_point() + theme_classic()
ggplot(data) + aes(x = Surface.Area.mm2, y = Volume.mm3, label = ids) + geom_point() + theme_classic() + geom_text()

ggplot(
  data %>% mutate(across(where(is.numeric), scale)) %>% pivot_longer(cols = data %>% select(where(is.numeric)) %>% colnames())
) + aes(x = name, y = value, label = ids) + geom_point() + 
  geom_text(data = . %>% filter(value >= 3 | value <= -3), position = position_nudge(x = 0.6)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90))

### Filtering of data for future tests

#What genotypes saw an odd number of strain-nema treatments? (note this filter should be before filtering by number of rhizobia strains)
data %>% ungroup %>% group_by(Genotype, Rhizo, Nema) %>% summarize(n = n()) %>% group_by(Genotype, .add = TRUE) %>% summarize(Nema_rhizo_combos = n()) %>% filter(Nema_rhizo_combos < 2)

step1_filter_data <- data %>%
  filter(!(Genotype == "HM041" & Rhizo == "H") & 
           !(Genotype == "HM044" & Rhizo == "WSM1022") & 
           !(Genotype == "HM049" & Rhizo == "WSM1022") & 
           !(Genotype == "HM065" & Rhizo == "I") & 
           !(Genotype == "HM148" & Rhizo == "H") & 
           !(Genotype == "HM167" & Rhizo == "I") & 
           !(Genotype == "HM167" & Rhizo == "WSM1022"))

#What genotypes saw only one rhizobia strain?
step1_filter_data %>% group_by(Genotype, Rhizo) %>% summarize() %>% group_by(Genotype, .add = TRUE) %>% summarize(Rhizo_strains = n()) %>% filter((Rhizo_strains < 2))

my_data <- step1_filter_data %>%
  filter(!Genotype %in% c("HM044", "HM168", "HM170")) %>%
  filter(!ids %in% c("id500", "id213"))


# Double checking filtering by applying checks for filtering again
my_data %>% group_by(Genotype, Rhizo, Nema) %>% summarize() %>% group_by(Genotype, .add = TRUE) %>% summarize(n = n()) %>% filter(n != 2)
my_data %>% group_by(Genotype, Rhizo) %>% summarize() %>% group_by(Genotype, .add = TRUE) %>% summarize(n = n()) %>% filter(n < 2)



### Setting subclusters for testing interaction terms
full_design <- my_data %>% select(Genotype, Rhizo) %>% mutate(value = 1) %>% distinct() %>% complete(Genotype, Rhizo)
full_design_2 <- my_data %>% group_by(Genotype, Rhizo) %>% summarize(samples = n()) %>% as.data.frame() %>% mutate(value = 1) %>% distinct() %>% complete(Genotype, Rhizo)
full_design$value <- replace(full_design$value, is.na(full_design$value), 0)
ggplot(full_design) + geom_tile(aes(x = Genotype, y = Rhizo, fill = value))
tibble <- full_design %>% pivot_wider(names_from = Genotype, values_from = value, values_fill = list(value = 0))
mat <- as.matrix(tibble[, -1])
colnames(mat) <- colnames(tibble)[2:ncol(tibble)]
rownames(mat) <- tibble$Rhizo


dist_matrix <- dist(mat, method = "binary")
hclust_result <- hclust(dist_matrix)

dist_matrix_cols <- dist(t(mat), method = "binary")
hclust_result_cols <- hclust(dist_matrix_cols)

Geno_order <- unique(full_design$Genotype)[hclust_result_cols$order]
Rhizo_order <- unique(full_design$Rhizo)[hclust_result$order]

ggplot(full_design) + geom_tile(aes(x = Genotype, y = Rhizo, fill = value)) +
  scale_x_discrete(limits = Geno_order) + scale_y_discrete(limits = Rhizo_order) + 
  theme(axis.text.x = element_text(angle = 90)) + theme(legend.position = "none")

ggplot(full_design) + geom_tile(aes(x = Genotype, y = Rhizo, fill = value, alpha = value)) +
  scale_x_discrete(limits = Geno_order) + scale_y_discrete(limits = Rhizo_order) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(data = full_design_2 %>% group_by(Rhizo) %>% summarize(n = sum(value, na.rm = TRUE)),
            aes(y = Rhizo, label = n), x = 21) +
  geom_text(data = full_design_2 %>% group_by(Genotype) %>% summarize(n = sum(value, na.rm = TRUE)),
            aes(x = Genotype, label = n), y = 11) +
  coord_cartesian(clip = "off") +
  scale_fill_gradient(low = "white", high = "black") +
  scale_alpha_identity() +
  theme_classic() + theme(legend.position = "blank") +
  theme(plot.margin = margin(t = 20, r = 20, l = 5, b = 5, "points")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_vline(xintercept = seq(0.5,20.5,1)) +
  geom_hline(yintercept = seq(0.5,10.5,1))

subcluster_1 <- my_data %>% filter(Rhizo %in% c("C", "D") & 
                                     Genotype %in% c("HM049", "HM138"))
subcluster_2 <- my_data %>% filter(Rhizo %in% c("F", "E") & 
                                     Genotype %in% c("HM012", "HM052"))
subcluster_3 <- my_data %>% filter(Rhizo %in% c("J", "I") & 
                                     Genotype %in% c("HM011", "HM071"))
subcluster_4 <- my_data %>% filter(Rhizo %in% c("E", "D") & 
                                     Genotype %in% c("HM012", "HM138"))
subcluster_5 <- my_data %>% filter(Rhizo %in% c("J", "USDA1021") & 
                                     Genotype %in% c("HM084", "HM071"))
manuscript_subcluster <- my_data %>% filter(Rhizo %in% c("WSM1022", "USDA1021") &
                                              Genotype %in% c("HM145"))

colnames(my_data %>% select(is.numeric))
### Correlations
corrplot(cor(
  my_data %>% select(is.numeric) %>%
    select(-c("Old...Galls", "Old...Nod.Count", "X.2",
              "Galls_manuscript", "Really.Old...Galls..pre.stain.", 
              "Block.Row", "Nods", "Nods.Removed", "Galls_eggless",
              "EggSacs", "Nodules")) %>% drop_na()
), method="color",
type="upper", order="hclust", 
addCoef.col = "black", # Add coefficient of correlation
tl.col="black", tl.srt=45, #Text label color and rotation
sig.level = 0.01, insig = "blank",number.digits = 2, number.cex = 0.6, 
# hide correlation coefficient on the principal diagonal
diag=FALSE
)


seed <- abs(round(rnorm(1)*1000, 0))