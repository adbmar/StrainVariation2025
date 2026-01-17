#Mac's PCA coordinates from his work on neutral variation in HapMap Genotypes
# data found at HapMap_pca.tsv

PCA_dir <- file.path(dir_main,"GenotypeSelection")
PCA <- read.csv(file.path(PCA_dir, "HapMap_pca.tsv"), sep = "\t", header = FALSE)
colnames(PCA) <- c("Genotype", "x", "y")
PCA <- PCA %>% dplyr::select(Genotype, x, y)

#Corlett's susceptibility data
# data found at GenotypeMeans_NodulesAndGalls.csv
gs_data <- read.csv(file.path(PCA_dir, "GenotypeMeans_NodulesAndGalls.csv"))
colnames(gs_data)[1] <- "Genotype"
gs_data <- gs_data %>% dplyr::select(Genotype, nod.R, Galls)

gs_data <- left_join(gs_data, PCA, by = "Genotype")
head(gs_data)

library(ggrepel)


gs_plot <- ggplot(gs_data) + aes(x = x, y = y, fill = nod.R, color = Galls, label = Genotype) +
  geom_text_repel(vjust = 0, nudge_y = 0.01, color = "black") +
  geom_point(shape = 21, stroke = 2, size = 3) +
  theme_classic() +
  scale_color_gradient(name = "Nematode susceptibility",
    low = "yellow", high = "red") +
  scale_fill_gradient(name = "Nodulation",
                       low = "cyan", high = "darkblue") +
  ylim(c(-0.25, 0.12)) + xlim(c(-0.1, 0.2))
gs_plot


setwd(dir_out)
ggsave(filename = "GenotypeSelection.png", plot = gs_plot, device = "png",
       width = 2160, height = 1440, units = "px", scale = 1.75)

gs_plot2 <- ggplot(gs_data %>% filter(Genotype %in% (my_data %>% pull(Genotype) %>% unique()))) + aes(x = x, y = y, fill = nod.R, color = Galls, label = Genotype) +
  geom_text_repel(vjust = 0, nudge_y = 0.01, color = "black") +
  geom_point(shape = 21, stroke = 2, size = 3) +
  theme_classic() +
  scale_color_gradient(name = "Nematode susceptibility",
                       low = "yellow", high = "red") +
  scale_fill_gradient(name = "Nodulation",
                      low = "cyan", high = "darkblue") +
  ylim(c(-0.25, 0.12)) + xlim(c(-0.1, 0.2))
gs_plot2

setwd(dir_out)
ggsave(filename = "GenotypeSelection2.png", plot = gs_plot2, device = "png",
       width = 2160, height = 1440, units = "px", scale = 1.75)


