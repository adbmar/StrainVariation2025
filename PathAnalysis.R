library(boot)
library(piecewiseSEM)

options(contrasts = c("contr.sum", "contr.poly"))


### Setting up  data
resistance_data$Genotype <- factor(resistance_data$Genotype)
resistance_data$Rhizo     <- factor(resistance_data$Rhizo)

resistance_data_path <- resistance_data %>% filter(!is.na(Total.Nodules))
resistance_data_path <- resistance_data_path %>% filter(!is.na(Below.ground.biomass))
resistance_data_path$Below.ground.biomass <- scale(resistance_data_path$Below.ground.biomass)
resistance_data_path <- resistance_data_path %>% dplyr::select(c(Genotype, Rhizo, Total.Nodules, Galls, Below.ground.biomass))
resistance_data_path$Genotype <- factor(resistance_data_path$Genotype)
resistance_data_path$Rhizo <- factor(resistance_data_path$Rhizo)
resistance_data_path$Below.ground.biomass <- scale(as.numeric(resistance_data_path$Below.ground.biomass))
resistance_data_path$Total.Nodules <- as.numeric(resistance_data_path$Total.Nodules)
resistance_data_path$Galls <- as.numeric(resistance_data_path$Galls)

##############################
####### PATH ANALYSIS  #######
##############################
####### Rhizo to Galls #######
##############################

#### Statistical Tests ######

m.model <- lm(Below.ground.biomass ~ Rhizo + Genotype, data = resistance_data_path,
              contrasts = list(Rhizo = "contr.sum", Genotype = "contr.sum"))
summary(m.model)
# plotResiduals(simulateResiduals(m.model))

y.model <- glm.nb(Galls ~ Rhizo + Genotype + Below.ground.biomass, data = resistance_data_path,
                  contrasts = list(Rhizo = "contr.sum", Genotype = "contr.sum"))
summary(y.model)
# plotResiduals(simulateResiduals(y.model))

sem_model_all <- psem(m.model, y.model)
summary(sem_model_all, standardize = "scale")

statistical_tests <- summary(sem_model_all, standardize = "scale")
statistical_tests <- statistical_tests$coefficients %>% as.data.frame()


#### Coefficients ######

resistance_data_path_gaus <- resistance_data_path
resistance_data_path_gaus$Galls <- scale(resistance_data_path_gaus$Galls)

m.model <- lm(Below.ground.biomass ~ Rhizo + Genotype, data = resistance_data_path_gaus,
              contrasts = list(Rhizo = "contr.sum", Genotype = "contr.sum"))
summary(m.model)
# plotResiduals(simulateResiduals(m.model))

y.model <- lm(Galls ~ Rhizo + Genotype + Below.ground.biomass, data = resistance_data_path_gaus,
              contrasts = list(Rhizo = "contr.sum", Genotype = "contr.sum"))
summary(y.model)
# plotResiduals(simulateResiduals(y.model))

sem_model_all <- psem(m.model, y.model)
summary(sem_model_all, standardize = "scale")

coefficients <- summary(sem_model_all, standardize = "scale")
coefficients <- coefficients$coefficients %>% as.data.frame()

path_results <- full_join(
  coefficients %>% dplyr::select(Response, Predictor, Estimate),
  statistical_tests %>% dplyr::select(Response, Predictor, P.Value),
  by = c("Response", "Predictor")) %>%
  mutate(sig = if_else(P.Value < 0.05, TRUE, FALSE))

ggplot(path_results) +
  geom_segment(data = . %>% filter(Predictor == "Rhizo" & Response == "Galls"),
               x = 0.1, xend = 0.85, y = 0, yend = 0,
               aes(color = sig),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"))) + 
  geom_text(data = . %>% filter(Predictor == "Rhizo" & Response == "Galls"),
            x = 0.5, y = 0.1,
            aes(label = paste0("β=",Estimate,"\np=",P.Value))) +
  geom_segment(data = . %>% filter(Predictor == "Rhizo" & Response == "Below.ground.biomass"),
               x = 0, xend = 0.25, y = 0.15, yend = 0.85,
               aes(color = sig),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"))) + 
  geom_text(data = . %>% filter(Predictor == "Rhizo" & Response == "Below.ground.biomass"),
            x = 0.15, y = 0.7,
            aes(label = paste0("β=",Estimate,"\np=",P.Value))) +
  geom_segment(data = . %>% filter(Predictor == "Below.ground.biomass" & Response == "Galls"),
               x = 0.75, xend = 1, yend = 0.15, y = 0.85,
               aes(color = sig),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"))) + 
  geom_text(data = . %>% filter(Predictor == "Below.ground.biomass" & Response == "Galls"),
            x = 0.85, y = 0.7,
            aes(label = paste0("β=",Estimate,"\np=",P.Value))) +
  geom_text(x = 0, y = 0, label = "Rhizo") +
  geom_text(x = 1, y = 0, label = "Galls") +
  geom_text(x = 0.5, y = 1, label = "Root Biomass") +
  theme_void() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgray")) +
  theme(panel.spacing = unit(2, "cm"))


path_results2 <- rbind(
  rbind(path_results %>% filter(Predictor == "Rhizo = MAG154") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG154"),
  rbind(path_results %>% filter(Predictor == "Rhizo = USDA1021") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "USDA1021"),
  rbind(path_results %>% filter(Predictor == "Rhizo = KH48e") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "KH48e"),
  rbind(path_results %>% filter(Predictor == "Rhizo = WSM1022") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "WSM1022"),
  rbind(path_results %>% filter(Predictor == "Rhizo = MAG177") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG177"),
  rbind(path_results %>% filter(Predictor == "Rhizo = MAG533") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG533"),
  rbind(path_results %>% filter(Predictor == "Rhizo = MAG282") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG282"),
  rbind(path_results %>% filter(Predictor == "Rhizo = MAG358") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG358"),
  rbind(path_results %>% filter(Predictor == "Rhizo = MAG540") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG540"),
  rbind(path_results %>% filter(Predictor == "Rhizo = KH35c") %>% mutate(Predictor = "Rhizo"),
        path_results %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "KH35c"))

prop_med <- path_results2 %>% 
  mutate(path = paste0(Predictor, "_to_", Response)) %>%
  dplyr::select(Rhizo, path, Estimate) %>%
  pivot_wider(names_from = path, values_from = Estimate) %>%
  mutate(PropMed = (as.numeric(Rhizo_to_Below.ground.biomass) * as.numeric(Below.ground.biomass_to_Galls)) / (as.numeric(Rhizo_to_Below.ground.biomass) * as.numeric(Below.ground.biomass_to_Galls) + as.numeric(Rhizo_to_Galls)))

ggplot(path_results2) +
  geom_segment(data = . %>% filter(Predictor == "Rhizo" & Response == "Galls"),
               x = 0.1, xend = 0.85, y = 0, yend = 0,
               aes(color = sig),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"))) + 
  geom_text(data = . %>% filter(Predictor == "Rhizo" & Response == "Galls"),
            x = 0.5, y = 0.1,
            aes(label = paste0("β=",Estimate,"\np=",P.Value))) +
  geom_segment(data = . %>% filter(Predictor == "Rhizo" & Response == "Below.ground.biomass"),
               x = 0, xend = 0.25, y = 0.15, yend = 0.85,
               aes(color = sig),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"))) + 
  geom_text(data = . %>% filter(Predictor == "Rhizo" & Response == "Below.ground.biomass"),
            x = 0.15, y = 0.7,
            aes(label = paste0("β=",Estimate,"\np=",P.Value))) +
  geom_segment(data = . %>% filter(Predictor == "Below.ground.biomass" & Response == "Galls"),
               x = 0.75, xend = 1, yend = 0.15, y = 0.85,
               aes(color = sig),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"))) + 
  geom_text(data = . %>% filter(Predictor == "Below.ground.biomass" & Response == "Galls"),
            x = 0.85, y = 0.7,
            aes(label = paste0("β=",Estimate,"\np=",P.Value))) +
  geom_text(x = 0, y = 0, label = "Rhizo") +
  geom_text(x = 1, y = 0, label = "Galls") +
  geom_text(x = 0.5, y = 1, label = "Root Biomass") +
  theme_void() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgray")) +
  facet_wrap(~Rhizo) +
  theme(panel.spacing = unit(2, "cm"))



boot_fun <- function(data, indices, verbose = FALSE){
  # Resample data
  d <- resistance_data_path_gaus[indices, ]
  if(verbose == TRUE){
    print(i)
    i <<- i + 1}
  
  boot_m.model <- lm(Below.ground.biomass ~ Total.Nodules + Rhizo + Genotype, data = d,
                     contrasts = list(Rhizo = "contr.sum", Genotype = "contr.sum"))
  boot_y.model <- lm(Galls ~ Total.Nodules + Rhizo + Genotype + Below.ground.biomass, data = d,
                     contrasts = list(Rhizo = "contr.sum", Genotype = "contr.sum"))
  boot_sem_model_all <- psem(boot_m.model, boot_y.model)
  boot_coefs_df <- summary(boot_sem_model_all,
                           standardize = "scale", 
                           test.type = "III",
                           intercepts = TRUE)$coefficients %>%
    dplyr::select(Response, Predictor, Estimate)
  
  boot_coefs_df2 <- rbind(
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = MAG154") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG154"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = USDA1021") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "USDA1021"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = KH48e") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "KH48e"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = WSM1022") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "WSM1022"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = MAG177") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG177"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = MAG533") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG533"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = MAG282") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG282"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = MAG358") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG358"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = MAG540") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "MAG540"),
    rbind(boot_coefs_df %>% filter(Predictor == "Rhizo = KH35c") %>% mutate(Predictor = "Rhizo"),
          boot_coefs_df %>% filter(Predictor == "Below.ground.biomass")) %>% mutate(Rhizo = "KH35c"))
  
  boot_out_df <- boot_coefs_df2 %>% 
    mutate(path = paste0(Predictor, "_to_", Response)) %>%
    dplyr::select(Rhizo, path, Estimate) %>%
    pivot_wider(names_from = path, values_from = Estimate) %>%
    mutate(IndirectPath = as.numeric(Rhizo_to_Below.ground.biomass) * as.numeric(Below.ground.biomass_to_Galls),
           DirectPath = as.numeric(Rhizo_to_Galls),
           TotalPath = as.numeric(Rhizo_to_Below.ground.biomass) * as.numeric(Below.ground.biomass_to_Galls) + as.numeric(Rhizo_to_Galls),
           PropMed = as.numeric(Rhizo_to_Below.ground.biomass) * as.numeric(Below.ground.biomass_to_Galls) / (as.numeric(Rhizo_to_Below.ground.biomass) * as.numeric(Below.ground.biomass_to_Galls) + as.numeric(Rhizo_to_Galls))
    ) %>% dplyr::select(Rhizo, Rhizo_to_Below.ground.biomass, Below.ground.biomass_to_Galls, Rhizo_to_Galls,
                        DirectPath, IndirectPath, TotalPath, PropMed)
  
  boot_out <- c(as.numeric(boot_out_df$Rhizo_to_Below.ground.biomass),
                as.numeric(boot_out_df$Below.ground.biomass_to_Galls),
                as.numeric(boot_out_df$Rhizo_to_Galls),
                as.numeric(boot_out_df$IndirectPath),
                as.numeric(boot_out_df$DirectPath),
                as.numeric(boot_out_df$TotalPath),
                as.numeric(boot_out_df$PropMed))
  
  names(boot_out) <- c(boot_out_df %>% rowwise() %>% mutate(get = paste0(Rhizo, ".RhizoToBiomass")) %>% pull(get),
                       boot_out_df %>% rowwise() %>% mutate(get = paste0(Rhizo, ".BiomassToGalls")) %>% pull(get),
                       boot_out_df %>% rowwise() %>% mutate(get = paste0(Rhizo, ".RhizoToGalls")) %>% pull(get),
                       boot_out_df %>% rowwise() %>% mutate(get = paste0(Rhizo, ".IndirectPath")) %>% pull(get),
                       boot_out_df %>% rowwise() %>% mutate(get = paste0(Rhizo, ".DirectPath")) %>% pull(get),
                       boot_out_df %>% rowwise() %>% mutate(get = paste0(Rhizo, ".TotalPath")) %>% pull(get),
                       boot_out_df %>% rowwise() %>% mutate(get = paste0(Rhizo, ".PropMed")) %>% pull(get))
  
  if(verbose == TRUE){print(boot_out)}
  return(boot_out)
}

boot_file <- file.path(dir_out, "bootstrapping_pathanalysis4_R2000.rds")

if(file.exists(boot_file)) {
  message("Cached path analysis boostrap file found. Loading...")
  boot_res <- readRDS(boot_file)
} else {
  message("No cached bootstrap file found for path analysis. Running bootstrap (this may take time)...")
  set.seed(123)
  i <- 1
  # boot_res <- boot(resistance_data_path, boot_fun, R = 100, verbose = TRUE)
  # boot_res <- boot(resistance_data_path, boot_fun, R = 2000, verbose = TRUE)
  boot_res <- boot(resistance_data_path, boot_fun, R = 10000, verbose = TRUE)
  
  
  saveRDS(boot_res, boot_file)
  message("Bootstrap complete. Results saved to file at the following location:")
  print(boot_file)
}

# View summary
boot_res


boot_df <- as.data.frame(boot_res$t)
colnames(boot_df) <- names(boot_res$t0)
ncol(boot_df)
nrow(boot_df)

# Compute mean, sd, and 95% CI for each strain
summary_df_long <- boot_df %>%
  summarise(across(everything(),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE),
                        lwr = ~quantile(., 0.025),
                        upr = ~quantile(., 0.975))
  )
  ) %>%
  tidyr::pivot_longer(everything(),
                      names_to = c("Rhizo", ".value"),
                      names_sep = c("_")) %>%
  rowwise() %>%
  mutate(sig = ifelse((lwr > 0 & upr > 0)|(lwr < 0 & upr < 0), TRUE, FALSE)) %>%
  mutate(Path = str_split(Rhizo, "\\.")[[1]][2],
         Rhizo = str_split(Rhizo, "\\.")[[1]][1])

summary_df_long %>% filter(sig) %>% arrange(Rhizo)

