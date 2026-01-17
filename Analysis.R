options(contrasts = c("contr.sum", "contr.poly"))
library(MASS)
library(MuMIn)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggpattern)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(colorspace)
library(piecewiseSEM)
library(grid)


#Definitions for visuzliaing
effect_order <- rev(c("Genotype", "Rhizo", "Nema", "Genotype:Rhizo", "Genotype.Rhizo", "Nema:Genotype", "Genotype.1", "Nema:Rhizo", "Rhizo.1", "Nema:Genotype:Rhizo", "Genotype.Rhizo.1", "Residuals"))
color_scheme <- c("Nema" = "firebrick2",
                  "Genotype" = "green4",
                  "Rhizo" = "dodgerblue2",
                  "Genotype:Rhizo" = "lightseagreen",
                  "Genotype.Rhizo" = "lightseagreen",
                  "Nema:Genotype" = "green4",
                  "Genotype.1" = "green4",
                  "Nema:Rhizo" = "dodgerblue2",
                  "Rhizo.1" = "dodgerblue2",
                  "Nema:Genotype:Rhizo" = "lightseagreen",
                  "Genotype.Rhizo.1" = "lightseagreen",
                  "Black" = "black")

geno_background <- desaturate(lighten(color_scheme["Genotype"], 0.85), 0.8)
rhizo_background <- desaturate(lighten(color_scheme["Rhizo"], 0.85), 0.8)
genorhizo_background <- desaturate(lighten(color_scheme["Genotype:Rhizo"], 0.85), 0.8)


rhizo_color_scheme_1 <- c(brewer.pal(name = "Blues", n = 6))
rhizo_color_scheme_1[1] <- "#202475"
rhizo_color_scheme_1[1] <- "#203575"
rhizo_color_scheme_2 <- desaturate(rhizo_color_scheme_1, amount = -0.15)
rhizo_color_scheme <- c(rhizo_color_scheme_1, rhizo_color_scheme_2)
rhizo_color_scheme <- rhizo_color_scheme[rhizo_color_scheme != "NA"]
names(rhizo_color_scheme) <- my_data %>% pull(Rhizo) %>% unique()
rhizo_color_scheme[2] <- "royalblue"
rhizo_color_scheme[3] <- "darkblue"
rhizo_color_scheme[8] <- "royalblue3"

geno_color_scheme_1 <- c(brewer.pal(name = "Greens", n = 7))
geno_color_scheme_2 <- desaturate(geno_color_scheme_1, amount = 0.3)
geno_color_scheme_3 <- desaturate(geno_color_scheme_1, amount = 0.6)
geno_color_scheme <- c(geno_color_scheme_1, geno_color_scheme_2, geno_color_scheme_3)
names(geno_color_scheme) <- data %>% pull(Genotype) %>% unique()
geno_color_scheme <- geno_color_scheme[!is.na(names(geno_color_scheme))]
geno_color_scheme[1] <- "chartreuse4"
geno_color_scheme[8] <- "darkgreen"
geno_color_scheme[9] <- "green4"
geno_color_scheme[15] <- "chartreuse3"
geno_color_scheme[16] <- "green3"

blank <- ggplot() + theme_void() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))

percent_modifier <- 0.01

#Functions
extract_variances_glmm <- function(model){
  results <- data.frame(
    "Effect" = character(0),
    "Variance" = numeric(0)
  )
  var_corr <- VarCorr(model)$cond
  for (effect in names(var_corr)) {
    results <- plyr::rbind.fill(results, data.frame(
      "Effect" = effect,
      "Variance" = var_corr[[effect]][1])
    )
  }
  results <- as.data.frame(results)
  rownames(results) <- results$Effect
  
  fam <- model$modelInfo$family$family
  # Predicted means on the response scale
  mu <- predict(model, type = "response")
  # Dispersion / size parameter
  disp <- sigma(model)
  if (fam == "gaussian") {
    resid_var <- disp^2
  } else if (fam == "lognormal") {
    # Fitted values on the log scale
    mu_log <- predict(model)
    sigma_log <- disp
    resid_var <- mean((exp(sigma_log^2) - 1) * exp(2 * mu_log + sigma_log^2))
  } else if (fam == "nbinom1") {
    # Variance = mu * (1 + alpha)
    alpha <- disp
    resid_var <- mean(mu * (1 + alpha))
  } else if (fam == "nbinom2") {
    # Variance = mu + mu^2 / theta
    theta <- disp
    resid_var <- mean(mu + mu^2 / theta)
  } else if (fam == "poisson") {
    # Variance = mu
    resid_var <- mean(mu)
  } else {
    stop("Family not supported for residual variance calculation.")
  }
  residual_df <- as.data.frame(t(c("Effect" = "Residuals", "Variance" = as.numeric(resid_var))))
  results <- plyr::rbind.fill(results, residual_df)
  results <- as.data.frame(results)
  results$Variance <- as.numeric(results$Variance)
  rownames(results) <- results$Effect
  return(results)
}

extract_blups_glmm <- function(model, effect){
  model_ranef <- ranef(model)
  results <- data.frame(
    "Effect" = character(0),
    "Level" = character(0),
    "BLUP" = numeric(),
    "condVar" = numeric(0)
  )
  for (effect in names(model_ranef$cond)){
    condVars <- attr(model_ranef$cond[[effect]], "condVar")
    new_results <- model_ranef$cond[[effect]] %>%
      as.data.frame() %>%
      rename("BLUP" = "(Intercept)") %>%
      rownames_to_column("Level") %>%
      mutate(Effect = effect,
             condVar = NA)
    new_results$condVar <- c(condVars)
    results <- rbind(results, new_results)
  }
  results <- results %>% rowwise() %>%
    mutate(Genotype = case_when(Effect == "Genotype" ~ Level,
                                Effect == "Genotype:Rhizo" ~ str_split(Level, ":")[[1]][1],
                                Effect == "Genotype.Rhizo" ~ str_split(Level, ".")[[1]][1],
                                Effect == "Nema:Genotype" ~ str_split(Level, ":")[[1]][2],
                                Effect == "Genotype.1" ~ str_split(Level, ".")[[1]][1],
                                Effect == "Nema:Genotype:Rhizo" ~ str_split(Level, ":")[[1]][2],
                                Effect == "Genotype.Rhizo.1" ~ str_split(Level, ":")[[1]][1],
                                .default = NA),
           Rhizo = case_when(Effect == "Rhizo" ~ Level,
                             Effect == "Genotype:Rhizo" ~ str_split(Level, ":")[[1]][2],
                             Effect == "Genotype.Rhizo" ~ str_split(Level, ".")[[1]][2],
                             Effect == "Nema:Rhizo" ~ str_split(Level, ":")[[1]][2],
                             Effect == "Rhizo.1" ~ str_split(Level, ".")[[1]][1],
                             Effect == "Nema:Genotype:Rhizo" ~ str_split(Level, ":")[[1]][3],
                             Effect == "Genotype.Rhizo.1" ~ str_split(Level, ".")[[1]][2],
                             .default = NA),
           Nema = case_when(Effect == "Nema" ~ Level,
                            Effect == "Nema:Genotype" ~ str_split(Level, ":")[[1]][1],
                            Effect == "Nema:Rhizo" ~ str_split(Level, ":")[[1]][1],
                            Effect == "Nema:Genotype:Rhizo" ~ str_split(Level, ":")[[1]][1],
                            .default = NA),
           Nema = ifelse(Nema == "+", "Infected", "Uninfected")
    )
  results <- as.data.frame(results)
  return(results)
}


halfway_vector <- function(vector){
  result <- numeric(length(vector))
  result[1] <- vector[1]/2
  for (i in 2:length(vector)){
    result[i] <- (vector[i-1] + vector[i])/2}
  return(result)
}

results_handling <- function(table){
  result_table <- table %>% filter(Effect != "Residuals")
  result_table$Effect <- factor(
    result_table$Effect,
    levels = effect_order)
  result_table <- result_table[rev(rownames(result_table)), ]
  result_table <- result_table %>% 
    mutate(lrt_pvalue = as.numeric(lrt_pvalue),
           VariancePercent = Variance/sum(Variance),
           mod_Variance = Variance + sum(Variance)*percent_modifier,
           mod_VariancePercent = mod_Variance/sum(mod_Variance), 
           Nema = ifelse(Effect %in% c("Nema:Genotype", "Nema:Rhizo", "Nema:Genotype:Rhizo",
                                       "Rhizo.1", "Genotype.1", "Genotype.Rhizo.1"),
                         "Infected", 
                         "Uninfected"),
           sig = paste("p =",round(lrt_pvalue, 3)),
           sig = ifelse(sig == "p = 0", "p < 0.001", sig)) %>%
    mutate(labels = case_when(Effect == "Genotype" ~ "Host genotype",
                              Effect == "Rhizo" ~ "Rhizobia strain",
                              Effect == "Genotype:Rhizo" ~ "Host genotype by rhizobia strain",
                              Effect == "Genotype.Rhizo" ~ "Host genotype by rhizobia strain",
                              Effect == "Nema:Genotype" ~ "Infection by host genotype",
                              Effect == "Genotype.1" ~ "Infection severity by host genotype",
                              Effect == "Nema:Rhizo" ~ "Infection by rhizobia strain",
                              Effect == "Rhizo.1" ~ "Infection severity by rhizobia strain",
                              Effect == "Nema:Genotype:Rhizo" ~ "Infection by host geotype by rhizobia strain",
                              Effect == "Genotype.Rhizo.1" ~ "Infection severity by host geotype by rhizobia strain",
                              .default = "error"),
           labels2 = paste0(labels, " (", scales::percent(VariancePercent, accuracy = 1),", ", sig, ")")) %>%
    slice(match(rev(effect_order), Effect)) %>%
    mutate(cumsum = cumsum(mod_VariancePercent),
           midway = halfway_vector(cumsum)) %>%
    as.data.frame()
  result_table <- plyr::rbind.fill(result_table, table %>% filter(Effect == "Residuals")) %>% 
    mutate(VariancePercent_withResids = Variance/sum(Variance),
           mod_Variance_withResids = Variance + sum(Variance)*percent_modifier,
           mod_VariancePercent_withResids = mod_Variance_withResids/sum(mod_Variance_withResids),
           labels2_withResids = paste0(labels, "(", scales::percent(VariancePercent_withResids, accuracy = 1)," ", sig, ")")) %>%
    slice(match(rev(effect_order), Effect)) %>%
    mutate(cumsum_withResids = cumsum(mod_VariancePercent_withResids),
           midway_withResids = halfway_vector(cumsum_withResids)) %>%
    as.data.frame()
  result_table$Effect <- factor(
    result_table$Effect,
    levels = effect_order)
  return(result_table)
}

add_in_subs <- function(table){
  return_table <- table %>%
    mutate(
      sub1 = ifelse(Genotype %in% (subcluster_1 %>% pull(Genotype)) &
                      Rhizo %in% (subcluster_1 %>% pull(Rhizo)),
                    TRUE, FALSE),
      sub2 = ifelse(Genotype %in% (subcluster_2 %>% pull(Genotype)) &
                      Rhizo %in% (subcluster_2 %>% pull(Rhizo)),
                    TRUE, FALSE),
      sub3 = ifelse(Genotype %in% (subcluster_3 %>% pull(Genotype)) &
                      Rhizo %in% (subcluster_3 %>% pull(Rhizo)),
                    TRUE, FALSE),
      sub4 = ifelse(Genotype %in% (subcluster_4 %>% pull(Genotype)) &
                      Rhizo %in% (subcluster_4 %>% pull(Rhizo)),
                    TRUE, FALSE),
      sub5 = ifelse(Genotype %in% (subcluster_5 %>% pull(Genotype)) &
                      Rhizo %in% (subcluster_5 %>% pull(Rhizo)),
                    TRUE, FALSE)) %>%
    pivot_longer(cols = c(sub1, sub2, sub3, sub4, sub5)) %>%
    rename("subcluster" = "name") %>%
    mutate(subcluster = case_when(subcluster == "sub1" ~ "Subcluster 1",
                                  subcluster == "sub2" ~ "Subcluster 2",
                                  subcluster == "sub3" ~ "Subcluster 3",
                                  subcluster == "sub4" ~ "Subcluster 4", 
                                  subcluster == "sub5" ~ "Subcluster 5",
                                  .default = "No subcluster")) %>%
    group_by(across(-c(subcluster, value))) %>%
    # group_by(Level, BLUP, Effect, condVar, Genotype, Rhizo, Nema) %>%
    reframe(subcluster = if (any(value)) subcluster[value] else NA) %>%
    unnest(subcluster) %>%
    as.data.frame()
  return(return_table)
}

# Setting up datasets
data_set_pre_filter <- data %>%
  filter(!is.na(Block.Position1)) %>%
  filter(!is.na(Block.Shelf)) %>%
  filter(!is.na(Rack))

virulence_data <- data_set_pre_filter %>% 
  filter(!is.na(Above.ground.biomass)) %>%
  filter(!is.na(Below.ground.biomass))

resistance_data <- data_set_pre_filter %>% 
  filter(Nema == "+") %>%
  filter(!is.na(Galls)) %>%
  filter(!is.na(Below.ground.biomass)) %>%
  filter(!is.na(Genotype)) %>%
  filter(!is.na(Rhizo))

resistance_data$Rhizo <- as.factor(resistance_data$Rhizo)
resistance_data$Genotype <- as.factor(resistance_data$Genotype)

tolerance_data <- data_set_pre_filter %>%
  filter(Nema == "+") %>%
  filter(!is.na(Galls)) %>%
  filter(!is.na(Below.ground.biomass)) %>%
  filter(!is.na(Above.ground.biomass))

mutualist_data <- data_set_pre_filter %>%
  filter(!is.na(Total.Nodules)) %>%
  filter(!is.na(Below.ground.biomass))

#######################
### ggplot formulas ###
#######################

VCA_plot <- theme_classic +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        legend.position = "none",
        plot.title.position = "plot") +
  scale_x_continuous(position = "top", labels = scales::percent) +
  scale_fill_manual(values = color_scheme) +
  scale_pattern_manual(values = c("Infected" = "stripe", "Uninfected" = "none")) +
  coord_cartesian(clip = "off", expand = 0)


#################
### Virulence ###
#################
{
  virulence_full_model <- glmmTMB(data = virulence_data,
                                  formula = Above.ground.biomass ~ Block.Position1 + Block.Shelf + Block.Side +
                                    Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo) + (1|Nema:Genotype:Rhizo),
                                  family = lognormal)
  # plot(testResiduals(virulence_full_model))
  plot(simulateResiduals(virulence_full_model))
  r.squaredGLMM(virulence_full_model)
  
  virulence_results <- extract_variances_glmm(virulence_full_model)
  virulence_results <- virulence_results %>% mutate(lrt_pvalue = "error")
  
  virulence_second_order_model <- virulence_full_model
  virulence_first_order_model <- virulence_full_model
  
  virulence_second_order_model <- update(virulence_second_order_model, .~.-(1|Nema:Genotype:Rhizo))
  virulence_first_order_model <- update(virulence_first_order_model, .~.-(1|Genotype:Rhizo)-(1|Nema:Rhizo)-(1|Nema:Genotype)-(1|Nema:Genotype:Rhizo))
  
  virulence_results["Genotype", "lrt_pvalue"] <- 
    anova(virulence_first_order_model,
          update(virulence_first_order_model, .~.-(1|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  virulence_results["Rhizo", "lrt_pvalue"] <- 
    anova(virulence_first_order_model,
          update(virulence_first_order_model, .~.-(1|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  virulence_results["Genotype:Rhizo", "lrt_pvalue"] <- 
    anova(virulence_second_order_model,
          update(virulence_second_order_model, .~.-(1|Genotype:Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  virulence_results["Nema:Genotype", "lrt_pvalue"] <- 
    anova(virulence_second_order_model,
          update(virulence_second_order_model, .~.-(1|Nema:Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  virulence_results["Nema:Rhizo", "lrt_pvalue"] <- 
    anova(virulence_second_order_model,
          update(virulence_second_order_model, .~.-(1|Nema:Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  virulence_results["Nema:Genotype:Rhizo", "lrt_pvalue"] <- 
    anova(virulence_full_model,
          virulence_second_order_model,
          type = "LRT")$`Pr(>Chisq)`[2]
  
  
  virulence_results_final <- results_handling(virulence_results)%>% rowwise() %>%
    mutate(labels2 = case_when(Effect == "Nema:Genotype" ~ paste("B.", labels2),
                               Effect == "Nema:Rhizo" ~ paste("C.", labels2),
                               Effect == "Nema:Genotype:Rhizo" ~ paste("D.", labels2),
                               .default = labels2),
           label_border = case_when(Effect == "Nema:Genotype" ~ Effect,
                                    Effect == "Nema:Rhizo" ~ Effect,
                                    Effect == "Nema:Genotype:Rhizo" ~ Effect,
                                    .default = "Black")) %>%
    as.data.frame()
  
  {
    virulence_results_final$labels_x <- c(0.2, 0.25, 0.325, 0.455, 0.8675, 0.95, NA)
    virulence_results_final$labels_y <- c(0.98,0.95,0.92,0.8875,0.8875,0.45, NA)
    
    
    virulence_BLUPs <- extract_blups_glmm(virulence_full_model) %>% filter(Effect %in% c("Nema:Rhizo", "Nema:Genotype", "Nema:Genotype:Rhizo"))
    virulence_BLUP_maxy <- ceiling(max(virulence_BLUPs$BLUP + virulence_BLUPs$condVar)*100)/100
    virulence_BLUP_miny <- floor(min(virulence_BLUPs$BLUP - virulence_BLUPs$condVar)*100)/100
    
    virulence_VCA <- ggplot(virulence_results_final %>% filter(Effect != "Residuals")) +
      aes(color = Effect, fill = Effect,
          x = mod_VariancePercent,
          pattern = Nema) +
      geom_col_pattern(position = "stack", aes(y = NA),
                       color = "white",
                       linewidth = 1.3,
                       pattern_fill = color_scheme["Nema"],
                       pattern_color = NA) +
      geom_label(data = . %>%
                   #Setting the labels to go in the bar chart vs. beneath
                   filter(Effect %in% c("Genotype", "Genotype:Rhizo")) %>%
                   mutate(labels2 = str_replace(labels2, " \\(", "\\\n\\(")),
                 inherit.aes = FALSE, hjust = 0.5,
                 position = "identity",
                 aes(label = labels2, y = 1, x = midway)) + 
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x.bottom = element_blank(),
            legend.position = "none",
            plot.title.position = "plot",
            plot.margin = margin(t=0.2,l=0.4,b=0,r=1, "cm")) +
      xlab("A. Variance explained: host shoot biomass") +
      scale_x_continuous(position = "top", labels = scales::percent) +
      scale_fill_manual(values = color_scheme) +
      scale_pattern_manual(values = c("Infected" = "stripe", "Uninfected" = "none")) +
      coord_cartesian(clip = "off", expand = 0)
    
    virulence_VCA_transition <- ggplot(virulence_results_final %>% filter(Effect != "Residuals") %>% 
                                         #Setting the labels to beneath the bar chart vs. within
                                         filter(!Effect %in% c("Genotype", "Genotype:Rhizo"))
                                       ) +
      annotate("rect", xmin = 0.1, xmax = 0.49, ymin = 0.48, ymax = 0.8875,
               fill = geno_background, color = color_scheme["Genotype"], linewidth = 2) +
      annotate("rect", xmin = 0.51, xmax = 0.9, ymin = 0.48, ymax = 0.8875,
               fill = rhizo_background, color = color_scheme["Rhizo"], linewidth = 2) +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.01, ymax = 0.45,
               fill = genorhizo_background, color = color_scheme["Genotype:Rhizo"], linewidth = 2) +
      geom_segment(inherit.aes = FALSE,
                   aes(x = midway,
                       xend = labels_x,
                       y = 1, 
                       yend = labels_y,
                       color = Effect)) +
      aes(y = labels_y,
          x = labels_x) +
      geom_label(aes(label = labels2, color = label_border),
                 label.r = unit(0.5, "lines"),
                 hjust = 0.925, 
                 fill = "white") +
      scale_fill_manual(values = color_scheme, guide = "none") +
      scale_color_manual(values = color_scheme, guide = "none") +
      xlim(c(0,1)) + ylim(c(0,1)) +
      coord_cartesian(expand = 0, clip = "off") +
      theme_void() +
      theme(plot.margin = margin(t=0,l=0.4,b=0,r=0.4, "cm"),
            title = element_blank(),
            plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2))
    
    virulence_geno <- ggplot(virulence_BLUPs %>% filter(Effect == "Nema:Genotype") )+
      aes(x = Nema,
          group = Genotype,
          color = Genotype,
          y = BLUP,
          ymax = BLUP+condVar,
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(virulence_BLUP_miny, virulence_BLUP_maxy)) +
      ylab("Host Shoot Biomass BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = geno_color_scheme, guide = "none") +
      theme_classic() +
      theme(plot.background = element_rect(fill = NA, color = NA),
            plot.title.position = "plot", 
            plot.margin = margin(t=0.25,l=0.5,b=0,r=0.5, "cm"))
    
    virulence_rhizo <- ggplot(virulence_BLUPs %>% filter(Effect == "Nema:Rhizo") )+
      aes(x = Nema,
          group = Rhizo,
          color = Rhizo,
          y = BLUP,
          ymax = BLUP+condVar,
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(virulence_BLUP_miny, virulence_BLUP_maxy)) +
      ylab("Host Shoot Biomass BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = rhizo_color_scheme, guide = "none") +
      theme_classic() +
      theme(plot.background = element_rect(fill = NA, color = NA),
            plot.title.position = "plot",
            plot.margin = margin(t=0.25,l=0.5,b=0,r=0.5, "cm"))
    
    virulence_genorhizo <- ggplot(virulence_BLUPs %>% 
                                    filter(Effect == "Nema:Genotype:Rhizo") %>%
                                    add_in_subs() %>%
                                    filter(!is.na(subcluster))) +
      aes(x = Nema, color = Rhizo, group = interaction(Genotype, Rhizo), shape = Genotype,
          y = BLUP, ymax = BLUP+condVar, ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      facet_wrap(~subcluster, nrow = 1) +
      ylim(c(virulence_BLUP_miny, virulence_BLUP_maxy)) +
      ylab("Host Shoot\nBiomass BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = rhizo_color_scheme) +
      scale_shape_manual(values = c(2,1,3,4,16,17,8)) +
      theme_classic() +
      coord_cartesian(clip = "off") + 
      theme(legend.position = "bottom",
            legend.title.position = "top",
            legend.background = element_rect(fill = NA, color = NA),
            legend.key = element_rect(fill = NA, color = NA),
            plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(t=1.25,l=0.15,b=0.05,r=0.15, "cm"))
    
    virulence_inset_plot <- ggplotGrob(
      plot_grid(ncol = 1, nrow = 2, rel_heights = c(3,4),
                plot_grid(nrow = 1, ncol = 4, rel_widths = c(1,4,4,1),
                          blank,
                          virulence_geno,
                          virulence_rhizo,
                          blank),
                virulence_genorhizo
      )
    )
    
    virulence_model_annotation <- ggplot() + theme_void() +
      theme(plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2)) +
      annotation_custom(xmin = 0, xmax = 1, ymin = 0.1, ymax = 0.9,
        grob = textGrob(x = 0, just = "left",
                        gp = gpar(fontsize = 8),
                        paste(
"Variance component analysis model
     Response variable: Above ground biomass (host fitness)
     Fixed effects: Blocking variables, Nematode infection status
     Random variables: Host genotype, rhizobia strain, host-genotype-by-rhizobia-strain,
                       infection-by-host-genotype, infection-by-rhizoba-strain, infection-by-host-genotype-by-rhizobia-strain\n
     Conditional coefficient of correlation:", round(r.squaredGLMM(virulence_full_model)[2], 4),"
     Marginal coefficient of correlation:", round(r.squaredGLMM(virulence_full_model)[1], 4)
    )))
    
    virulence_title <- ggdraw() + 
      draw_label("Parasite virulence", fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7),
            plot.background = element_rect(fill = "white", color = "white"))
    
    virulence_plot <- plot_grid(ncol = 1, nrow = 4, align = "v", rel_heights = c(1,5,18,4),
                                virulence_title,
                                virulence_VCA,
                                virulence_VCA_transition +
                                  annotation_custom(
                                    virulence_inset_plot,
                                    xmin = 0, xmax = 1,
                                    ymin = 0, ymax = 0.875),
                                virulence_model_annotation
    )
    # virulence_plot
    ggsave(plot = virulence_plot, file = file.path(dir_out, "VirulencePlot.png"),
           width = 11880, height = 12060, unit = "px", dpi = 1200, scale = 0.9)
    }
}



#################
### Resistance ###
#################
{
  resistance_full_model <- glmmTMB(data = resistance_data,
                                   formula = Galls ~ scale(Below.ground.biomass) + Block.Position1 + Block.Shelf + Block.Side +
                                     (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo),
                                   family = nbinom1)
  # plot(testResiduals(resistance_full_model))
  # plot(simulateResiduals(resistance_full_model))
  
  resistance_results <- extract_variances_glmm(resistance_full_model)
  resistance_results <- resistance_results %>% mutate(lrt_pvalue = "error")
  
  resistance_second_order_model <- resistance_full_model
  resistance_first_order_model <- resistance_full_model
  
  resistance_second_order_model <- update(resistance_second_order_model, .~.-(1|Nema:Genotype:Rhizo))
  resistance_first_order_model <- update(resistance_first_order_model, .~.-(1|Genotype:Rhizo)-(1|Nema:Rhizo)-(1|Nema:Genotype)-(1|Nema:Genotype:Rhizo))
  
  resistance_results["Genotype", "lrt_pvalue"] <- 
    anova(resistance_full_model,
          update(resistance_first_order_model, .~.-(1|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  resistance_results["Rhizo", "lrt_pvalue"] <- 
    anova(resistance_full_model,
          update(resistance_first_order_model, .~.-(1|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  resistance_results["Genotype:Rhizo", "lrt_pvalue"] <- 
    anova(resistance_full_model,
          update(resistance_second_order_model, .~.-(1|Genotype:Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  
  
  resistance_results_final <- results_handling(resistance_results) %>% rowwise() %>%
    mutate(labels2 = case_when(Effect == "Genotype" ~ paste("B.", labels2),
                               Effect == "Rhizo" ~ paste("C.", labels2),
                               Effect == "Genotype:Rhizo" ~ paste("D.", labels2),
                               .default = labels2),
           label_border = Effect) %>% as.data.frame()
  
  {
    resistance_results_final$labels_x <- c(0.465, 0.85, 0.965, NA)
    resistance_results_final$labels_y <- c(0.95,0.95,0.45, NA)
    
    
    resistance_BLUPs <- extract_blups_glmm(resistance_full_model)
    resistance_BLUP_maxy <- ceiling(max(resistance_BLUPs$BLUP + resistance_BLUPs$condVar)*100)/100
    resistance_BLUP_miny <- floor(min(resistance_BLUPs$BLUP + resistance_BLUPs$condVar)*100)/100
    
    resistance_VCA <- ggplot(resistance_results_final %>% filter(Effect != "Residuals")) +
      aes(color = Effect, fill = Effect,
          x = mod_VariancePercent, 
          pattern = Nema) +
      geom_col_pattern(position = "stack", aes(y = NA),
                       color = "white",
                       linewidth = 1.3,
                       pattern_fill = color_scheme["Nema"],
                       pattern_color = NA) +
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.text.x = element_text(hjust = 0),
            legend.position = "none",
            plot.margin = margin(t = 0.2, l=0.4, b=0, r=1, "cm")) +
      xlab("A. Variance explained: host shoot biomass") +
      scale_x_continuous(position = "top", labels = scales::percent) +
      scale_fill_manual(values = color_scheme) +
      scale_pattern_manual(values = c("Infected" = "stripe", "Uninfected" = "none")) +
      coord_cartesian(clip = "off", expand = 0)
    
    resistance_VCA_transition <- ggplot(resistance_results_final %>% filter(Effect != "Residuals")) +
      annotate("rect", xmin = 0.125, xmax = 0.49, ymin = 0.48, ymax = 0.95,
               fill = geno_background, color = color_scheme["Genotype"], linewidth = 2) +
      annotate("rect", xmin = 0.51, xmax = 0.875, ymin = 0.48, ymax = 0.95,
               fill = rhizo_background, color = color_scheme["Rhizo"], linewidth = 2) +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.01, ymax = 0.45,
               fill = genorhizo_background, color = color_scheme["Genotype:Rhizo"], linewidth = 2) +
      geom_segment(inherit.aes = FALSE,
                   aes(x = midway, 
                       xend = labels_x,
                       y = 1,
                       yend = labels_y,
                       color = Effect)) +
      aes(y = labels_y,
          x = labels_x) +
      geom_label(aes(label = labels2, color = label_border),
                 label.r = unit(0.5, "lines"),
                 hjust = 0.925,
                 fill = "white") +
      scale_fill_manual(values = color_scheme, guide = "none") +
      scale_color_manual(values = color_scheme, guide = "none") +
      xlim(c(0,1)) + ylim(c(0,1)) +
      coord_cartesian(expand = 0, clip = "off") +
      theme_void() +
      theme(plot.margin = margin(t=0,l=0.4,b=0,r=0.4, "cm"),
            title = element_blank(),
            plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2))
    
    resistance_geno <- ggplot(resistance_BLUPs %>% filter(Effect == "Genotype") )+
      aes(x = Genotype, 
          color = Genotype,
          y = BLUP,
          ymax = BLUP+condVar, 
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange() +
      theme_classic() +
      ylim(c(resistance_BLUP_miny, resistance_BLUP_maxy)) +
      ylab("Number of Galls BLUP") +
      xlab("Host genotype") +
      scale_color_manual(values = geno_color_scheme, guide = "none") +
      theme(plot.background = element_rect(fill = NA, color = NA),
            axis.text.x = element_text(angle = 90),
            plot.margin = margin(t=0.75, l=0.5, b=0, r=0.5, "cm"))
    
    resistance_rhizo <- ggplot(resistance_BLUPs %>% filter(Effect == "Rhizo") )+
      aes(x = Rhizo, 
          color = Rhizo,
          y = BLUP, 
          ymax = BLUP+condVar, 
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(resistance_BLUP_miny, resistance_BLUP_maxy)) +
      ylab("Number of Galls BLUP") +
      xlab("Rhizobia strain") +
      scale_color_manual(values = rhizo_color_scheme, guide = "none") +
      theme_classic() +
      theme(plot.background = element_rect(fill = NA, color = NA),
            axis.text.x = element_text(angle = 90),
            plot.margin = margin(t=0.75, l=0.5, b=0, r=0.5, "cm"))
    
    resistance_genorhizo <- ggplot(resistance_BLUPs %>%
                                     filter(Effect == "Genotype:Rhizo") %>%
                                     add_in_subs() %>%
                                     filter(!is.na(subcluster))) +
      aes(x = Genotype, color = Rhizo, group = Rhizo,
          y = BLUP, ymax = BLUP+condVar, ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(resistance_BLUP_miny, resistance_BLUP_maxy)) +
      facet_wrap(~subcluster, nrow = 1, scales = "free_x") +
      ylab("Number of Galls BLUP") +
      xlab("Host genotype") +
      scale_color_manual(values = rhizo_color_scheme) +
      theme_classic() +
      theme(
        legend.title.position = "left",
        plot.background = element_rect(fill = NA, color = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = "bottom",
        legend.justification.bottom = c(0,0),
        plot.margin = margin(t=1.25, l=0.25, b=0.25, r=0.25, "cm")) +
      guides(color = guide_legend(nrow = 1))
    
    resistance_inset_plot <- ggplotGrob(
      plot_grid(ncol = 1, nrow = 2,
                plot_grid(ncol = 4, nrow = 1, rel_widths = c(1,3,3,1), align = "h",
                          blank,
                          resistance_geno,
                          resistance_rhizo,
                          blank),
                resistance_genorhizo
      )
    )
    
    
    resistance_model_annotation <- ggplot() + theme_void() +
      theme(plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2)) +
      annotation_custom(xmin = 0, xmax = 1, ymin = 0.1, ymax = 0.9,
                        grob = textGrob(x = 0, just = "left",
                                        gp = gpar(fontsize = 8),
                                        paste(
"Variance component analysis model
     Response variable: Gall count
     Fixed effects: Blocking variables, Root biomass,
     Random variables: Host genotype, rhizobia strain, host-genotype-by-rhizobia-strain\n
     Conditional coefficient of correlation:", round(r.squaredGLMM(resistance_full_model)[6], 4),"
     Marginal coefficient of correlation:", round(r.squaredGLMM(resistance_full_model)[3], 4)
                                        )))
 
    resistance_title <- ggdraw() + 
      draw_label("Host resistance", fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7),
            plot.background = element_rect(fill = "white", color = "white"))
    
    resistance_plot <- plot_grid(ncol = 1, nrow = 4, align = "v", rel_heights = c(1,5,18,4), 
                                 resistance_title,
                                 resistance_VCA,
                                 resistance_VCA_transition +
                                   annotation_custom(
                                     resistance_inset_plot,
                                     xmin = 0, xmax = 1,
                                     ymin = 0, ymax = 0.975),
                                 resistance_model_annotation
    )
    # resistance_plot
    ggsave(plot = resistance_plot, file = file.path(dir_out, "ResistancePlot.png"),
           width = 11880, height = 12060, unit = "px", dpi = 1200, scale = 0.9)
           }
}


#################
### Mutualist ###
#################
{
  mutualist_full_model <- glmmTMB(data = mutualist_data,
                                  formula = Total.Nodules ~ scale(Below.ground.biomass) + Block.Position1 + Block.Shelf + Block.Side +
                                    Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo) + (1|Nema:Genotype:Rhizo),
                                  family = nbinom1)
  # plot(testResiduals(mutualist_full_model))
  # plot(simulateResiduals(mutualist_full_model))
  
  mutualist_results <- extract_variances_glmm(mutualist_full_model)
  mutualist_results <- mutualist_results %>% mutate(lrt_pvalue = "error")
  
  mutualist_first_order_model <- mutualist_full_model
  mutualist_second_order_model <- mutualist_full_model
  
  mutualist_second_order_model <- update(mutualist_second_order_model, .~.-(1|Nema:Genotype:Rhizo))
  mutualist_first_order_model <- update(mutualist_first_order_model, .~.-(1|Genotype:Rhizo)-(1|Nema:Rhizo)-(1|Nema:Genotype)-(1|Nema:Genotype:Rhizo))
  
  mutualist_results["Genotype", "lrt_pvalue"] <- 
    anova(mutualist_first_order_model,
          update(mutualist_first_order_model, .~.-(1|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  mutualist_results["Rhizo", "lrt_pvalue"] <- 
    anova(mutualist_first_order_model,
          update(mutualist_first_order_model, .~.-(1|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  mutualist_results["Genotype:Rhizo", "lrt_pvalue"] <- 
    anova(mutualist_second_order_model,
          update(mutualist_second_order_model, .~.-(1|Genotype:Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  mutualist_results["Nema:Genotype", "lrt_pvalue"] <- 
    anova(mutualist_second_order_model,
          update(mutualist_second_order_model, .~.-(1|Nema:Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  mutualist_results["Nema:Rhizo", "lrt_pvalue"] <- 
    anova(mutualist_second_order_model,
          update(mutualist_second_order_model, .~.-(1|Nema:Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  mutualist_results["Nema:Genotype:Rhizo", "lrt_pvalue"] <- 
    anova(mutualist_full_model,
          mutualist_second_order_model,
          type = "LRT")$`Pr(>Chisq)`[2]
  
  
  mutualist_results_final <- results_handling(mutualist_results) %>% rowwise() %>%
    mutate(labels2 = case_when(Effect == "Nema:Genotype" ~ paste("B.", labels2),
                               Effect == "Nema:Rhizo" ~ paste("C.", labels2),
                               Effect == "Nema:Genotype:Rhizo" ~ paste("D.", labels2),
                               .default = labels2),
           label_border = case_when(Effect == "Nema:Genotype" ~ Effect,
                                    Effect == "Nema:Rhizo" ~ Effect,
                                    Effect == "Nema:Genotype:Rhizo" ~ Effect,
                                    .default = "Black")) %>%
    as.data.frame()
  
  {
    mutualist_results_final$labels_x <- c(0.2, 0.25, 0.325, 0.45, 0.845, 0.95, NA)
    mutualist_results_final$labels_y <- c(0.98,0.95,0.92,0.8875,0.8875,0.45, NA)
    
    mutualist_BLUPs <- extract_blups_glmm(mutualist_full_model) %>% filter(Effect %in% c("Nema:Genotype", "Nema:Rhizo", "Nema:Genotype:Rhizo"))
    mutualist_BLUP_maxy <- ceiling(max(mutualist_BLUPs$BLUP + mutualist_BLUPs$condVar)*100)/100
    mutualist_BLUP_miny <- floor(min(mutualist_BLUPs$BLUP + mutualist_BLUPs$condVar)*100)/100
    
    mutualist_VCA <- ggplot(mutualist_results_final %>% filter(Effect != "Residuals")) +
      aes(color = Effect, fill = Effect, 
          x = mod_VariancePercent,
          pattern = Nema) +
      geom_col_pattern(position = "stack", aes(y = NA),
                       color = "white", 
                       linewidth = 1.3,
                       pattern_fill = color_scheme["Nema"],
                       pattern_color = NA) +
      geom_label(data = . %>%
                   #Setting the labels to go in the bar chart vs. beneath
                   filter(Effect %in% c("Genotype", "Genotype:Rhizo")) %>%
                   mutate(labels2 = str_replace(labels2, " \\(", "\\\n\\("),
                          labels2 = str_replace(labels2, "by ", "by\\\n")),
                 inherit.aes = FALSE, hjust = 0.5,
                 position = "identity",
                 aes(label = labels2, y = 1, x = midway)) + 
      xlab("A. Variance explained: number of nodules") +
      scale_x_continuous(position = "top", labels = scales::percent) +
      scale_fill_manual(values = color_scheme) +
      scale_pattern_manual(values = c("Infected" = "stripe", "Uninfected" = "none")) +
      theme(plot.margin = margin(t=0.75,l=0.4,b=0,r=1, "cm")) +
      coord_cartesian(clip = "off", expand = 0) +
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.text.x = element_text(hjust = 0),
            legend.position = "none",
            plot.margin = margin(t=0.75,l=0.4,b=0,r=0.75, "cm"))
    
    mutualist_VCA_transition <- ggplot(mutualist_results_final %>% filter(Effect != "Residuals") %>% 
                                         #Setting the labels to beneath the bar chart vs. within
                                         filter(!Effect %in% c("Genotype", "Genotype:Rhizo"))
    ) +
      annotate("rect", xmin = 0.125, xmax = 0.49, ymin = 0.48, ymax = 0.8875,
               fill = geno_background, color = color_scheme["Genotype"], linewidth = 2) +
      annotate("rect", xmin = 0.51, xmax = 0.875, ymin = 0.48, ymax = 0.8875,
               fill = rhizo_background, color = color_scheme["Rhizo"], linewidth = 2) +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.01, ymax = 0.45,
               fill = genorhizo_background, color = color_scheme["Genotype:Rhizo"], linewidth = 2) +
      geom_segment(inherit.aes = FALSE,
                   aes(x = midway, 
                       xend = labels_x,
                       y = 1,
                       yend = labels_y,
                       color = Effect)) +
      aes(y = labels_y, 
          x = labels_x) +
      geom_label(aes(label = labels2, color = label_border),
                 label.r = unit(0.5, "lines"),
                 hjust = 0.925,
                 fill = "white") +
      scale_fill_manual(values = color_scheme, guide = "none") +
      scale_color_manual(values = color_scheme, guide = "none") +
      xlim(c(0,1)) + ylim(c(0,1)) +
      coord_cartesian(expand = 0, clip = "off") +
      theme_void() +
      theme(plot.margin = margin(t=0, l=0.4, b=0, r=0.4, "cm"),
            title = element_blank(),
            plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2))
    
    mutualist_geno <- ggplot(mutualist_BLUPs %>% filter(Effect == "Nema:Genotype") )+
      aes(x = Nema,
          group = Genotype, 
          color = Genotype,
          y = BLUP, 
          ymax = BLUP+condVar,
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(mutualist_BLUP_miny, mutualist_BLUP_maxy)) +
      ylab("Number of Nodules BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = geno_color_scheme, guide = "none") +
      theme_classic() +
      theme(plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(t=0.25, l=0.5, b=0, r=0.5, "cm"))
    
    mutualist_rhizo <- ggplot(mutualist_BLUPs %>% filter(Effect == "Nema:Rhizo") )+
      aes(x = Nema, 
          group = Rhizo,
          color = Rhizo,
          y = BLUP,
          ymax = BLUP+condVar, 
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(mutualist_BLUP_miny, mutualist_BLUP_maxy)) +
      ylab("Number of Nodules BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = rhizo_color_scheme, guide = "none") +
      theme_classic() +
      theme(plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(t=0.25, l=0.5, b=0, r=1, "cm"))
    
    mutualist_genorhizo <- ggplot(mutualist_BLUPs %>% 
                                    filter(Effect == "Nema:Genotype:Rhizo") %>%
                                    add_in_subs() %>%
                                    filter(!is.na(subcluster))) +
      aes(x = Nema, color = Rhizo, group = interaction(Genotype, Rhizo), shape = Genotype,
          y = BLUP, ymax = BLUP+condVar, ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      facet_wrap(~subcluster, nrow = 1) +
      ylim(c(mutualist_BLUP_miny, mutualist_BLUP_maxy)) +
      ylab("Number of\nNodules BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = rhizo_color_scheme) +
      scale_shape_manual(values = c(2,1,3,4,16,17,8)) +
      theme_classic() +
      theme(legend.position = "bottom",
            legend.title.position = "top",
            legend.background = element_rect(fill = NA, color = NA),
            legend.key = element_rect(fill = NA, color = NA),
            plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(t=1.25,l=0.15,b=0.05,r=0.15, "cm"))
    
    
    mutualist_inset_plot <- ggplotGrob(
      plot_grid(ncol = 1, nrow = 2, rel_heights = c(3,4),
                plot_grid(nrow = 1, ncol = 4, rel_widths = c(1,3,3,1),
                          blank,
                          mutualist_geno,
                          mutualist_rhizo,
                          blank),
                mutualist_genorhizo
      )
    )
    
    mutualist_model_annotation <- ggplot() + theme_void() +
      theme(plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2)) +
      annotation_custom(xmin = 0, xmax = 1, ymin = 0.1, ymax = 0.9,
                        grob = textGrob(x = 0, just = "left",
                                        gp = gpar(fontsize = 8),
                                        paste(
"Variance component analysis model
     Response variable: Nodule counts
     Fixed effects: Blocking variables, Root biomass,
     Random variables: Host genotype, rhizobia strain, host-genotype-by-rhizobia-strain,
                       infection-by-host-genotype, infection-by-rhizoba-strain, infection-by-host-genotype-by-rhizobia-strain\n
     Conditional coefficient of correlation:", round(r.squaredGLMM(mutualist_full_model)[6], 4),"
     Marginal coefficient of correlation:", round(r.squaredGLMM(mutualist_full_model)[3], 4)
                                        )))
    
    mutualist_title <- ggdraw() + 
      draw_label("Mutualism robustness", fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7),
            plot.background = element_rect(fill = "white", color = "white"))    
    
    mutualist_plot <- plot_grid(ncol = 1, nrow = 4, align = "v", rel_heights = c(1,5,18,4), 
                                mutualist_title,
                                mutualist_VCA,
                                mutualist_VCA_transition +
                                  annotation_custom(
                                    mutualist_inset_plot,
                                    xmin = 0, xmax = 1,
                                    ymin = 0, ymax = 0.875),
                                mutualist_model_annotation
    ) + theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm"),
              plot.background = element_rect(fill = "white", color = "white"))
    # mutualist_plot
    ggsave(plot = mutualist_plot, file = file.path(dir_out, "MutualistPlot.png"),
           width = 11880, height = 12060, unit = "px", dpi = 1200, scale = 0.9)
    }
}



#################
### Tolerance ###
#################
{
  tolerance_data_mod <- tolerance_data %>% mutate(Below.ground.biomass = scale(Below.ground.biomass), Galls = Galls)
  tolerance_full_model <- glmmTMB(data = tolerance_data_mod,
                                  formula = log(Above.ground.biomass) ~ Below.ground.biomass + Block.Position1 + Block.Shelf + Block.Side +
                                    Galls + (1|Genotype) + (1|Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo),
                                  family = gaussian)
  r.squaredGLMM(tolerance_full_model)
  # plot(testResiduals(tolerance_full_model))
  plot(simulateResiduals(tolerance_full_model))
  
  tolerance_results <- extract_variances_glmm(tolerance_full_model)
  tolerance_results <- tolerance_results %>% mutate(lrt_pvalue = "error") %>% as.data.frame()
  
  tolerance_first_order_model <- tolerance_full_model
  tolerance_first_order_model <- update(tolerance_first_order_model, .~.-(0+Galls|Rhizo)-(0+Galls|Genotype))
  
  tolerance_results["Genotype", "lrt_pvalue"] <- 
    anova(tolerance_first_order_model,
          update(tolerance_first_order_model, .~.-(1|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  tolerance_results["Rhizo", "lrt_pvalue"] <- 
    anova(tolerance_first_order_model,
          update(tolerance_first_order_model, .~.-(1|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  tolerance_results["Genotype.1", "lrt_pvalue"] <- 
    anova(tolerance_first_order_model,
          update(tolerance_full_model, .~.-(0+Galls|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  tolerance_results["Rhizo.1", "lrt_pvalue"] <- 
    anova(tolerance_first_order_model,
          update(tolerance_full_model, .~.-(0+Galls|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  
  
  tolerance_results_final <- results_handling(tolerance_results) %>% rowwise() %>%
    mutate(labels2 = case_when(Effect == "Genotype.1" ~ paste("B.", labels2),
                               Effect == "Rhizo.1" ~ paste("C.", labels2),
                               .default = labels2),
           label_border = case_when(Effect %in% c("Genotype.1", "Rhizo.1") ~ Effect,
                                    .default = "Black")) %>%
    as.data.frame()
  tolerance_results_final
  
  tolerance_results_final$labels_x <- c(0.2, 0.25, 0.455, 0.935, NA)
  tolerance_results_final$labels_y <- c(0.98,0.97,0.8875, 0.8875, NA)
  
  tolerance_VCA <- ggplot(tolerance_results_final %>% filter(Effect != "Residuals")) +
    aes(color = Effect, fill = Effect,
        x = mod_VariancePercent,
        pattern = Nema) +
    geom_col_pattern(position = "stack", aes(y = NA),
                     color = "white",
                     linewidth = 1.3,
                     pattern_fill = color_scheme["Nema"],
                     pattern_color = NA) +
    geom_label(data = . %>%
                 #Setting the labels to go in the bar chart vs. beneath
                 filter(Effect %in% c("Genotype")) %>%
                 mutate(labels2 = str_replace(labels2, " \\(", "\\\n\\("),
                        labels2 = str_replace(labels2, "by ", "by\n")),
               inherit.aes = FALSE, hjust = 0.5,
               position = "identity",
               aes(label = labels2, y = 1, x = midway)) + 
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x.bottom = element_blank(),
          axis.text.x = element_text(hjust = 0),
          legend.position = "none",
          plot.margin = margin(t=0.2,l=0.4,b=0,r=1, "cm")) +
    xlab("A. Variance explained: host shoot biomass") +
    scale_x_continuous(position = "top", labels = scales::percent) +
    scale_fill_manual(values = color_scheme) +
    scale_pattern_manual(values = c("Infected" = "stripe", "Uninfected" = "none")) +
    theme(plot.margin = margin(t=0.2,l=0.4,b=0,r=0.4, "cm")) +
    coord_cartesian(clip = "off", expand = 0)
  
  tolerance_VCA_transition <- ggplot(tolerance_results_final %>% filter(Effect != "Residuals") %>% 
                                       #Setting the labels to beneath the bar chart vs. within
                                       filter(!Effect %in% c("Genotype"))
  ) +
    annotate("rect", xmin = 0.025, xmax = 0.49, ymin = 0.05, ymax = 0.8875,
             fill = geno_background, color = color_scheme["Genotype"], linewidth = 2) +
    annotate("rect", xmin = 0.51, xmax = 0.975, ymin = 0.05, ymax = 0.8875,
             fill = rhizo_background, color = color_scheme["Rhizo"], linewidth = 2) +
    geom_segment(inherit.aes = FALSE,
                 aes(x = midway,
                     xend = labels_x,
                     y = 1, 
                     yend = labels_y,
                     color = Effect)) +
    aes(y = labels_y,
        x = labels_x) +
    geom_label(aes(label = labels2, color = label_border),
               label.r = unit(0.5, "lines"),
               hjust = 0.925,
               fill = "white") +
    scale_fill_manual(values = color_scheme, guide = "none") +
    scale_color_manual(values = color_scheme, guide = "none") +
    xlim(c(0,1)) + ylim(c(0,1)) +
    coord_cartesian(expand = 0, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t=0,l=0.4,b=0,r=0.4, "cm"),
          title = element_blank(),
          plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
          panel.background = element_rect(fill = "white", color = "white", linewidth = 2))
  
  # tolerance_residuals_model <- update(tolerance_full_model,
  # .~.-(Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo)))
  tolerance_residuals_model <- update(tolerance_full_model, Galls~Below.ground.biomass, family = gaussian)
  
  resids <- residuals(tolerance_residuals_model, type = "response")
  length(resids)
  dim(tolerance_data)
  tolerance_residual_data <- cbind(tolerance_data, "Residuals" = resids) %>% as.data.frame()
  
  cor(cbind("Galls" = tolerance_residual_data$Galls, "GallDensity" = tolerance_residual_data$Residuals, "ShootBiomass" = tolerance_residual_data$Above.ground.biomass))
  
  # tolerance_vis_model <- glmmTMB(data = tolerance_residual_data,
  # formula = Residuals ~ Galls + Genotype + Rhizo + (1|Genotype:Rhizo) + Galls:Genotype + Galls:Rhizo,
  # family = gaussian)
  # tolerance_vis_model_GR <- glmmTMB(data = tolerance_residual_data %>% mutate(GR = interaction(Genotype, Rhizo)),
  # formula = Residuals ~ Galls + GR + Galls:GR,
  # family = gaussian)
  # tolerance_vis_model <- glmmTMB(data = tolerance_residual_data,
  #                                formula = Above.ground.biomass ~ Residuals + Genotype + Rhizo + (1|Genotype:Rhizo) + Residuals:Genotype + Residuals:Rhizo,
  #                                family = lognormal)
  # tolerance_vis_model_GR <- glmmTMB(data = tolerance_residual_data %>% mutate(GR = interaction(Genotype, Rhizo)),
  #                                   formula = Above.ground.biomass ~ Residuals + GR + Residuals:GR,
  #                                   family = lognormal)
  
  tolerance_vis_model <- update(tolerance_full_model,
                                data = tolerance_residual_data,
                                formula = . ~ . - Below.ground.biomass - Galls - (1|Genotype) - (1|Rhizo) -
                                  (0+Galls|Genotype) - (0+Galls|Rhizo) -
                                  Block.Position1 - Block.Shelf - Block.Side +
                                  Residuals + Genotype + Rhizo + Residuals:Genotype + Residuals:Rhizo)
  
  plot(simulateResiduals(tolerance_vis_model))
  
  residual_minmax <- plyr::rbind.fill(
    tolerance_residual_data %>% group_by(Genotype, Rhizo) %>% summarize(min_residuals = min(Residuals), max_residuals = max(Residuals)),
    tolerance_residual_data %>% group_by(Rhizo) %>% summarize(min_residuals = min(Residuals), max_residuals = max(Residuals)),
    tolerance_residual_data %>% group_by(Genotype) %>% summarize(min_residuals = min(Residuals), max_residuals = max(Residuals))
  ) %>% mutate(Effect = ifelse(is.na(Rhizo) | is.na(Genotype),"Main" ,"Interaction"))
  
  
  tolerance_viz_rhizo <- emmip(tolerance_vis_model,
                               ~Residuals:Rhizo,
                               var = "Galls", 
                               specs = c("Rhizo", "Residuals"),
                               type = "response",
                               CIs = TRUE,
                               at = list(Residuals = c(seq(-2,4,0.2))),
                               plotit = FALSE, 
                               rg.limit = 1000000) %>% as.data.frame()
  
  tolerance_viz_geno <- emmip(tolerance_vis_model,
                              ~Residuals:Genotype,
                              var = "Residuals", 
                              specs = c("Genotype", "Residuals"),
                              type = "response", 
                              at = list(Residuals = c(seq(-2,4,0.2))),
                              CIs = TRUE,
                              plotit = FALSE, 
                              rg.limit = 1000000) %>% as.data.frame()
  
  
  tolerance_viz_rhizo_limited <- tolerance_viz_rhizo %>%
    left_join(residual_minmax %>% 
                filter(Effect == "Main", !is.na(Rhizo)),
              by = "Rhizo") %>%
    filter(Residuals >= min_residuals, Residuals <= max_residuals)
  
  tolerance_viz_geno_limited <- tolerance_viz_geno %>%
    left_join(residual_minmax %>% filter(Effect == "Main", !is.na(Genotype)),
              by = "Genotype") %>%
    filter(Residuals >= min_residuals, Residuals <= max_residuals)
  
  max_list <- sort(rbind(tolerance_viz_rhizo_limited, tolerance_viz_geno_limited)$UCL)
  tolerance_max <- ceiling(100 * max_list[length(max_list)])/100
  # tolerance_max <- 0.5
  
  tolerance_min <- floor(100 *min(rbind(tolerance_viz_rhizo_limited, tolerance_viz_geno_limited)$LCL, na.rm = TRUE))/100
  
  tolerance_rhizo <- ggplot(tolerance_viz_rhizo_limited) +
    aes(x = Residuals,
        y = yvar,
        ymax = UCL,
        ymin = LCL,
        color = Rhizo, 
        fill = Rhizo) + 
    geom_ribbon(alpha = 0.05, color = NA) +
    geom_line() +
    ylim(c(tolerance_min, tolerance_max)) +
    scale_color_manual(values = rhizo_color_scheme) +
    scale_fill_manual(values = rhizo_color_scheme) +
    ylab("Host Shoot Biomass\nMarginal Means") +
    xlab("Residual number of galls") +
    theme_classic() +
    theme(plot.background = element_rect(fill = NA, color = NA),
          legend.position = "none",
          plot.margin = margin(t = 0.25, b = 1, l = 0.5, r = 0, unit = "cm"))
  
  
  tolerance_geno <- ggplot(tolerance_viz_geno_limited %>% rowwise() %>% mutate(UCL = ifelse(UCL > tolerance_max, tolerance_max, UCL))) +
    aes(x = Residuals, y = yvar,
        ymax = UCL, ymin = LCL,
        color = Genotype, fill = Genotype) + 
    geom_ribbon(alpha = 0.05, color = NA) +
    geom_line() +
    ylim(c(tolerance_min, tolerance_max)) +
    scale_color_manual(values = geno_color_scheme) +
    scale_fill_manual(values = geno_color_scheme) +
    ylab("Host Shoot Biomass\nMarginal Means") +
    xlab("Residual number of galls") +
    theme_classic() +
    coord_cartesian(clip = "off") + 
    theme(plot.background = element_rect(fill = NA, color = NA),
          legend.position = "none",
          plot.margin = margin(t = 0.25, b = 1, l = 0, r = 0.5, unit = "cm"))
  

  
  
  tolerance_inset_plot <- ggplotGrob(
              plot_grid(nrow = 1, ncol = 4, rel_widths = c(1,12,12,1),
                        blank,
                        tolerance_geno,
                        tolerance_rhizo,
                        blank)
    )
  
  
  tolerance_model_annotation <- ggplot() + theme_void() +
    theme(plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
          panel.background = element_rect(fill = "white", color = "white", linewidth = 2)) +
    annotation_custom(xmin = 0, xmax = 1, ymin = 0.1, ymax = 0.9,
                      grob = textGrob(x = 0, just = "left",
                                      gp = gpar(fontsize = 8),
                                      paste(
"Variance component analysis model
     Response variable: Shoot biomass
     Fixed effects: Blocking variables, Gall counts, Root biomass,
     Random variables: Host genotype, rhizobia strain, host-genotype-by-rhizobia-strain,
                       gall-count-by-host-genotype, gall-count-by-rhizoba-strain, gall-count-by-host-genotype-by-rhizobia-strain\n
     Conditional coefficient of correlation:", round(r.squaredGLMM(tolerance_full_model)[2], 4),"
     Marginal coefficient of correlation:", round(r.squaredGLMM(tolerance_full_model)[1], 4)
                                      )))
  
  tolerance_title <- ggdraw() + 
    draw_label("Host tolerance", fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7),
          plot.background = element_rect(fill = "white", color = "white"))
  
  
  tolerance_plot <- plot_grid(ncol = 1, nrow = 4, align = "v", rel_heights = c(1,5,12,4), 
                              tolerance_title,
                              tolerance_VCA,
                              tolerance_VCA_transition +
                                annotation_custom(
                                  tolerance_inset_plot,
                                  xmin = 0, xmax = 1,
                                  ymin = 0, ymax = 0.875),
                              tolerance_model_annotation
  ) + theme(plot.margin = margin(r=0.5, unit = "cm"),
            plot.background = element_rect(fill = "white", color = "white"))
  # tolerance_plot
  ggsave(plot = tolerance_plot, file = file.path(dir_out, "TolerancePlot.png"),
         width = 11880, height = 8040, unit = "px", dpi = 1200, scale = 1)
}

# tolerance_plot
# virulence_plot
# resistance_plot
# mutualist_plot

r.squaredGLMM(resistance_full_model)
r.squaredGLMM(virulence_full_model)
r.squaredGLMM(tolerance_full_model)
r.squaredGLMM(mutualist_full_model)







#################
### Resistance ###
#################
{
  nobgb_resistance_full_model <- glmmTMB(data = resistance_data,
                                         formula = Galls ~ Block.Position1 + Block.Shelf + Block.Side +
                                           (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo),
                                         family = nbinom1)
  # plot(testResiduals(nobgb_resistance_full_model))
  # plot(simulateResiduals(nobgb_resistance_full_model))
  
  nobgb_resistance_results <- extract_variances_glmm(nobgb_resistance_full_model)
  nobgb_resistance_results <- nobgb_resistance_results %>% mutate(lrt_pvalue = "error")
  
  nobgb_resistance_second_order_model <- update(nobgb_resistance_full_model, .~.-(1|Nema:Genotype:Rhizo))
  nobgb_resistance_first_order_model <- update(nobgb_resistance_second_order_model, .~.-(1|Genotype:Rhizo)-(1|Nema:Rhizo)-(1|Nema:Genotype))
  
  nobgb_resistance_results["Genotype", "lrt_pvalue"] <- 
    anova(nobgb_resistance_full_model,
          update(nobgb_resistance_first_order_model, .~.-(1|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_resistance_results["Rhizo", "lrt_pvalue"] <- 
    anova(nobgb_resistance_full_model,
          update(nobgb_resistance_first_order_model, .~.-(1|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_resistance_results["Genotype:Rhizo", "lrt_pvalue"] <- 
    anova(nobgb_resistance_full_model,
          update(nobgb_resistance_second_order_model, .~.-(1|Genotype:Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  
  
  nobgb_resistance_results_final <- results_handling(nobgb_resistance_results) %>% rowwise() %>%
    mutate(labels2 = case_when(Effect == "Genotype" ~ paste("B.", labels2),
                               Effect == "Rhizo" ~ paste("C.", labels2),
                               Effect == "Genotype:Rhizo" ~ paste("D.", labels2),
                               .default = labels2),
           label_border = Effect) %>% as.data.frame()
  
  {
    nobgb_resistance_results_final$labels_x <- c(0.465, 0.85, 0.965, NA)
    nobgb_resistance_results_final$labels_y <- c(0.975,0.975,0.45, NA)
    
    
    nobgb_resistance_BLUPs <- extract_blups_glmm(nobgb_resistance_full_model)
    nobgb_resistance_BLUP_maxy <- ceiling(max(nobgb_resistance_BLUPs$BLUP + nobgb_resistance_BLUPs$condVar)*100)/100
    nobgb_resistance_BLUP_miny <- floor(min(nobgb_resistance_BLUPs$BLUP + nobgb_resistance_BLUPs$condVar)*100)/100
    
    nobgb_resistance_VCA <- ggplot(nobgb_resistance_results_final %>% filter(Effect != "Residuals")) +
      aes(color = Effect, fill = Effect,
          x = mod_VariancePercent, 
          pattern = Nema) +
      geom_col_pattern(position = "stack", aes(y = NA),
                       color = "white",
                       linewidth = 1.3,
                       pattern_fill = color_scheme["Nema"],
                       pattern_color = NA) +
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.text.x = element_text(hjust = 0),
            legend.position = "none",
            plot.margin = margin(t = 0.2, l=0.4, b=0, r=1, "cm")) +
      xlab("A. Variance components: host shoot biomass") +
      scale_x_continuous(position = "top", labels = scales::percent) +
      scale_fill_manual(values = color_scheme) +
      scale_pattern_manual(values = c("Infected" = "stripe", "Uninfected" = "none")) +
      coord_cartesian(clip = "off", expand = 0)
    
    nobgb_resistance_VCA_transition <- ggplot(nobgb_resistance_results_final %>% filter(Effect != "Residuals")) +
      annotate("rect", xmin = 0.125, xmax = 0.49, ymin = 0.48, ymax = 0.975,
               fill = geno_background, color = color_scheme["Genotype"], linewidth = 2) +
      annotate("rect", xmin = 0.51, xmax = 0.875, ymin = 0.48, ymax = 0.975,
               fill = rhizo_background, color = color_scheme["Rhizo"], linewidth = 2) +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.01, ymax = 0.45,
               fill = genorhizo_background, color = color_scheme["Genotype:Rhizo"], linewidth = 2) +
      geom_segment(inherit.aes = FALSE,
                   aes(x = midway, 
                       xend = labels_x,
                       y = 1,
                       yend = labels_y,
                       color = Effect)) +
      aes(y = labels_y,
          x = labels_x) +
      geom_label(aes(label = labels2, color = label_border),
                 label.r = unit(0.5, "lines"),
                 hjust = 0.925,
                 fill = "white") +
      scale_fill_manual(values = color_scheme, guide = "none") +
      scale_color_manual(values = color_scheme, guide = "none") +
      xlim(c(0,1)) + ylim(c(0,1)) +
      coord_cartesian(expand = 0, clip = "off") +
      theme_void() +
      theme(plot.margin = margin(t=0,l=0.4,b=0,r=0.4, "cm"),
            title = element_blank(),
            plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2))
    
    nobgb_resistance_geno <- ggplot(nobgb_resistance_BLUPs %>% filter(Effect == "Genotype") )+
      aes(x = Genotype, 
          color = Genotype,
          y = BLUP,
          ymax = BLUP+condVar, 
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange() +
      theme_classic() +
      ylim(c(nobgb_resistance_BLUP_miny, nobgb_resistance_BLUP_maxy)) +
      ylab("Number of Galls BLUP") +
      xlab("Host genotype") +
      scale_color_manual(values = geno_color_scheme, guide = "none") +
      theme(plot.background = element_rect(fill = NA, color = NA),
            axis.text.x = element_text(angle = 90),
            plot.margin = margin(t=0.5, l=0.5, b=0, r=0.5, "cm"))
    
    nobgb_resistance_rhizo <- ggplot(nobgb_resistance_BLUPs %>% filter(Effect == "Rhizo") )+
      aes(x = Rhizo, 
          color = Rhizo,
          y = BLUP, 
          ymax = BLUP+condVar, 
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(nobgb_resistance_BLUP_miny, nobgb_resistance_BLUP_maxy)) +
      ylab("Number of Galls BLUP") +
      xlab("Rhizobia strain") +
      scale_color_manual(values = rhizo_color_scheme, guide = "none") +
      theme_classic() +
      theme(plot.background = element_rect(fill = NA, color = NA),
            axis.text.x = element_text(angle = 90),
            plot.margin = margin(t=0.5, l=0.5, b=0, r=0.5, "cm"))
    
    nobgb_resistance_genorhizo <- ggplot(nobgb_resistance_BLUPs %>%
                                           filter(Effect == "Genotype:Rhizo") %>%
                                           add_in_subs() %>%
                                           filter(!is.na(subcluster))) +
      aes(x = Genotype, color = Rhizo, group = Rhizo,
          y = BLUP, ymax = BLUP+condVar, ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(nobgb_resistance_BLUP_miny, nobgb_resistance_BLUP_maxy)) +
      facet_wrap(~subcluster, nrow = 1, scales = "free_x") +
      ylab("Number of Galls BLUP") +
      xlab("Host genotype") +
      scale_color_manual(values = rhizo_color_scheme) +
      theme_classic() +
      theme(
        legend.title.position = "left",
        plot.background = element_rect(fill = NA, color = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = "bottom",
        legend.justification.bottom = c(0,0),
        plot.margin = margin(t=1.25, l=0.25, b=0.25, r=0.25, "cm")) +
      guides(color = guide_legend(nrow = 1))
    
    nobgb_resistance_inset_plot <- ggplotGrob(
      plot_grid(ncol = 1, nrow = 2,
                plot_grid(ncol = 4, nrow = 1, rel_widths = c(1,3,3,1), align = "h",
                          blank,
                          nobgb_resistance_geno,
                          nobgb_resistance_rhizo,
                          blank),
                nobgb_resistance_genorhizo
      )
    )
    
    nobgb_resistance_model_annotation <- ggplot() + theme_void() +
      theme(plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2)) +
      annotation_custom(xmin = 0, xmax = 1, ymin = 0.1, ymax = 0.9,
                        grob = textGrob(x = 0, just = "left",
                                        gp = gpar(fontsize = 8),
                                        paste(
"Variance component analysis model
     Response variable: Gall counts
     Fixed effects: Blocking variables
     Random variables: Host genotype, rhizobia strain, host-genotype-by-rhizobia-strain\n
     Conditional coefficient of correlation:", round(r.squaredGLMM(nobgb_resistance_full_model)[6], 4),"
     Marginal coefficient of correlation:", round(r.squaredGLMM(nobgb_resistance_full_model)[3], 4)
                                        )))
    
    
    nobgb_resistance_plot <- plot_grid(ncol = 1, nrow = 4, align = "v", rel_heights = c(1,5,18,4), 
                                       resistance_title,
                                       nobgb_resistance_VCA,
                                       nobgb_resistance_VCA_transition +
                                         annotation_custom(
                                           nobgb_resistance_inset_plot,
                                           xmin = 0, xmax = 1,
                                           ymin = 0, ymax = 0.975),
                                       nobgb_resistance_model_annotation
    )
    # nobgb_resistance_plot
    ggsave(plot = nobgb_resistance_plot, file = file.path(dir_out, "ResistancePlot_nobgb.png"),
           width = 11880, height = 12060, unit = "px", dpi = 1200, scale = 0.9)
           }
}


#################
### Mutualist ###
#################
{
  nobgb_mutualist_full_model <- glmmTMB(data = mutualist_data,
                                        formula = Total.Nodules ~ Block.Position1 + Block.Shelf + Block.Side +
                                          Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo) + (1|Nema:Genotype:Rhizo),
                                        family = nbinom1)
  # plot(testResiduals(nobgb_mutualist_full_model))
  # plot(simulateResiduals(nobgb_mutualist_full_model))
  
  nobgb_mutualist_results <- extract_variances_glmm(nobgb_mutualist_full_model)
  nobgb_mutualist_results <- nobgb_mutualist_results %>% mutate(lrt_pvalue = "error")
  
  nobgb_mutualist_second_order_model <- update(nobgb_mutualist_full_model, .~.-(1|Nema:Genotype:Rhizo))
  nobgb_mutualist_first_order_model <- update(nobgb_mutualist_second_order_model, .~.-(1|Genotype:Rhizo)-(1|Nema:Rhizo)-(1|Nema:Genotype))
  
  nobgb_mutualist_results["Genotype", "lrt_pvalue"] <- 
    anova(nobgb_mutualist_first_order_model,
          update(nobgb_mutualist_first_order_model, .~.-(1|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_mutualist_results["Rhizo", "lrt_pvalue"] <- 
    anova(nobgb_mutualist_first_order_model,
          update(nobgb_mutualist_first_order_model, .~.-(1|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_mutualist_results["Genotype:Rhizo", "lrt_pvalue"] <- 
    anova(nobgb_mutualist_second_order_model,
          update(nobgb_mutualist_second_order_model, .~.-(1|Genotype:Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_mutualist_results["Nema:Genotype", "lrt_pvalue"] <- 
    anova(nobgb_mutualist_second_order_model,
          update(nobgb_mutualist_second_order_model, .~.-(1|Nema:Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_mutualist_results["Nema:Rhizo", "lrt_pvalue"] <- 
    anova(nobgb_mutualist_second_order_model,
          update(nobgb_mutualist_second_order_model, .~.-(1|Nema:Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_mutualist_results["Nema:Genotype:Rhizo", "lrt_pvalue"] <- 
    anova(nobgb_mutualist_full_model,
          nobgb_mutualist_second_order_model,
          type = "LRT")$`Pr(>Chisq)`[2]
  
  
  nobgb_mutualist_results_final <- results_handling(nobgb_mutualist_results) %>% rowwise() %>%
    mutate(labels2 = case_when(Effect == "Nema:Genotype" ~ paste("B.", labels2),
                               Effect == "Nema:Rhizo" ~ paste("C.", labels2),
                               Effect == "Nema:Genotype:Rhizo" ~ paste("D.", labels2),
                               .default = labels2),
           label_border = case_when(Effect == "Nema:Genotype" ~ Effect,
                                    Effect == "Nema:Rhizo" ~ Effect,
                                    Effect == "Nema:Genotype:Rhizo" ~ Effect,
                                    .default = "Black")) %>%
    as.data.frame()
  
  {
    nobgb_mutualist_results_final$labels_x <- c(0.2, 0.25, 0.325, 0.45, 0.845, 0.95, NA)
    nobgb_mutualist_results_final$labels_y <- c(0.98,0.95,0.92,0.8875,0.8875,0.45, NA)
    
    nobgb_mutualist_BLUPs <- extract_blups_glmm(nobgb_mutualist_full_model) %>% filter(Effect %in% c("Nema:Genotype", "Nema:Rhizo", "Nema:Genotype:Rhizo"))
    nobgb_mutualist_BLUP_maxy <- ceiling(max(nobgb_mutualist_BLUPs$BLUP + nobgb_mutualist_BLUPs$condVar)*100)/100
    nobgb_mutualist_BLUP_miny <- floor(min(nobgb_mutualist_BLUPs$BLUP + nobgb_mutualist_BLUPs$condVar)*100)/100
    
    nobgb_mutualist_VCA <- ggplot(nobgb_mutualist_results_final %>% filter(Effect != "Residuals")) +
      aes(color = Effect, fill = Effect, 
          x = mod_VariancePercent,
          pattern = Nema) +
      geom_col_pattern(position = "stack", aes(y = NA),
                       color = "white", 
                       linewidth = 1.3,
                       pattern_fill = color_scheme["Nema"],
                       pattern_color = NA) +
      xlab("A. Variance explained: number of nodules") +
      scale_x_continuous(position = "top", labels = scales::percent) +
      scale_fill_manual(values = color_scheme) +
      scale_pattern_manual(values = c("Infected" = "stripe", "Uninfected" = "none")) +
      theme(plot.margin = margin(t=0.2,l=0.4,b=0,r=1, "cm")) +
      coord_cartesian(clip = "off", expand = 0) +
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.text.x = element_text(hjust = 0),
            legend.position = "none",
            plot.margin = margin(t=0.2,l=0.4,b=0,r=0.75, "cm"))
    
    nobgb_mutualist_VCA_transition <- ggplot(nobgb_mutualist_results_final %>% filter(Effect != "Residuals")) +
      annotate("rect", xmin = 0.125, xmax = 0.49, ymin = 0.48, ymax = 0.8875,
               fill = geno_background, color = color_scheme["Genotype"], linewidth = 2) +
      annotate("rect", xmin = 0.51, xmax = 0.875, ymin = 0.48, ymax = 0.8875,
               fill = rhizo_background, color = color_scheme["Rhizo"], linewidth = 2) +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.01, ymax = 0.45,
               fill = genorhizo_background, color = color_scheme["Genotype:Rhizo"], linewidth = 2) +
      geom_segment(inherit.aes = FALSE,
                   aes(x = midway, 
                       xend = labels_x,
                       y = 1,
                       yend = labels_y,
                       color = Effect)) +
      aes(y = labels_y, 
          x = labels_x) +
      geom_label(aes(label = labels2, color = label_border),
                 label.r = unit(0.5, "lines"),
                 hjust = 0.925,
                 fill = "white") +
      scale_fill_manual(values = color_scheme, guide = "none") +
      scale_color_manual(values = color_scheme, guide = "none") +
      xlim(c(0,1)) + ylim(c(0,1)) +
      coord_cartesian(expand = 0, clip = "off") +
      theme_void() +
      theme(plot.margin = margin(t=0, l=0.4, b=0, r=0.4, "cm"),
            title = element_blank(),
            plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2))
    
    nobgb_mutualist_geno <- ggplot(nobgb_mutualist_BLUPs %>% filter(Effect == "Nema:Genotype") )+
      aes(x = Nema,
          group = Genotype, 
          color = Genotype,
          y = BLUP, 
          ymax = BLUP+condVar,
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(nobgb_mutualist_BLUP_miny, nobgb_mutualist_BLUP_maxy)) +
      ylab("Number of Nodules BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = geno_color_scheme, guide = "none") +
      theme_classic() +
      theme(plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(t=0.25, l=0.5, b=0, r=0.5, "cm"))
    
    nobgb_mutualist_rhizo <- ggplot(nobgb_mutualist_BLUPs %>% filter(Effect == "Nema:Rhizo") )+
      aes(x = Nema, 
          group = Rhizo,
          color = Rhizo,
          y = BLUP,
          ymax = BLUP+condVar, 
          ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      ylim(c(nobgb_mutualist_BLUP_miny, nobgb_mutualist_BLUP_maxy)) +
      ylab("Number of Nodules BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = rhizo_color_scheme, guide = "none") +
      theme_classic() +
      theme(plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(t=0.25, l=0.5, b=0, r=1, "cm"))
    
    nobgb_mutualist_genorhizo <- ggplot(nobgb_mutualist_BLUPs %>% 
                                          filter(Effect == "Nema:Genotype:Rhizo") %>%
                                          add_in_subs() %>%
                                          filter(!is.na(subcluster))) +
      aes(x = Nema, color = Rhizo, group = interaction(Genotype, Rhizo), shape = Genotype,
          y = BLUP, ymax = BLUP+condVar, ymin = BLUP-condVar) +
      geom_hline(yintercept = 0) +
      geom_pointrange(position = position_dodge(width = 0.25)) +
      geom_line(position = position_dodge(width = 0.25)) +
      facet_wrap(~subcluster, nrow = 1) +
      ylim(c(nobgb_mutualist_BLUP_miny, nobgb_mutualist_BLUP_maxy)) +
      ylab("Number of Nodules BLUP") +
      xlab("Nematode Infection Status") +
      scale_color_manual(values = rhizo_color_scheme) +
      scale_shape_manual(values = c(2,1,3,4,16,17,8)) +
      theme_classic() +
      theme(legend.position = "bottom",
            legend.title.position = "top",
            legend.background = element_rect(fill = NA, color = NA),
            legend.key = element_rect(fill = NA, color = NA),
            plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(t=1.25,l=0.15,b=0.05,r=0.15, "cm"))
    
    
    nobgb_mutualist_inset_plot <- ggplotGrob(
      plot_grid(ncol = 1, nrow = 2, rel_heights = c(3,4),
                plot_grid(nrow = 1, ncol = 4, rel_widths = c(1,3,3,1),
                          blank,
                          nobgb_mutualist_geno,
                          nobgb_mutualist_rhizo,
                          blank),
                nobgb_mutualist_genorhizo
      )
    )
    
    nobgb_mutualist_model_annotation <- ggplot() + theme_void() +
      theme(plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
            panel.background = element_rect(fill = "white", color = "white", linewidth = 2)) +
      annotation_custom(xmin = 0, xmax = 1, ymin = 0.1, ymax = 0.9,
                        grob = textGrob(x = 0, just = "left",
                                        gp = gpar(fontsize = 8),
                                        paste(
"Variance component analysis model
     Response variable: Nodule counts
     Fixed effects: Blocking variables,
     Random variables: Host genotype, rhizobia strain, host-genotype-by-rhizobia-strain,
                       infection-by-host-genotype, infection-by-rhizoba-strain, infection-by-host-genotype-by-rhizobia-strain\n
     Conditional coefficient of correlation:", round(r.squaredGLMM(nobgb_mutualist_full_model)[6], 4),"
     Marginal coefficient of correlation:", round(r.squaredGLMM(nobgb_mutualist_full_model)[3], 4)
                                        )))
    
    nobgb_mutualist_plot <- plot_grid(ncol = 1, nrow = 4, align = "v", rel_heights = c(1,5,18,4), 
                                      mutualist_title,
                                      nobgb_mutualist_VCA,
                                      nobgb_mutualist_VCA_transition +
                                        annotation_custom(
                                          nobgb_mutualist_inset_plot,
                                          xmin = 0, xmax = 1,
                                          ymin = 0, ymax = 0.875),
                                      nobgb_mutualist_model_annotation
    )
    # nobgb_mutualist_plot
    ggsave(plot = nobgb_mutualist_plot, file = file.path(dir_out, "MutualistPlot_nobgb.png"),
           width = 11880, height = 12060, unit = "px", dpi = 1200, scale = 0.9)
    }
}



#################
### Tolerance ###
#################
{
  nobgb_tolerance_data_mod <- tolerance_data %>% mutate(Below.ground.biomass = scale(Below.ground.biomass), Galls = Galls)
  nobgb_tolerance_full_model <- update(tolerance_full_model, .~.-Below.ground.biomass)
  # glmmTMB(data = nobgb_tolerance_data_mod,
  # formula = Above.ground.biomass ~ Galls + #Block.Position1 + Block.Shelf + Block.Side +
  # (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo),
  # family = lognormal)
  # plot(testResiduals(nobgb_tolerance_full_model))
  plot(simulateResiduals(nobgb_tolerance_full_model))
  
  nobgb_tolerance_results <- extract_variances_glmm(nobgb_tolerance_full_model)
  nobgb_tolerance_results <- nobgb_tolerance_results %>% mutate(lrt_pvalue = "error") %>% as.data.frame()
  
  nobgb_tolerance_first_order_model <- update(nobgb_tolerance_full_model, .~.-(0+Galls|Rhizo)-(0+Galls|Genotype))
  
  nobgb_tolerance_results["Genotype", "lrt_pvalue"] <- 
    anova(nobgb_tolerance_first_order_model,
          update(nobgb_tolerance_first_order_model, .~.-(1|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_tolerance_results["Rhizo", "lrt_pvalue"] <- 
    anova(nobgb_tolerance_first_order_model,
          update(nobgb_tolerance_first_order_model, .~.-(1|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_tolerance_results["Genotype.1", "lrt_pvalue"] <- 
    anova(nobgb_tolerance_full_model,
          update(nobgb_tolerance_full_model, .~.-(0+Galls|Genotype)),
          type = "LRT")$`Pr(>Chisq)`[2]
  nobgb_tolerance_results["Rhizo.1", "lrt_pvalue"] <- 
    anova(nobgb_tolerance_full_model,
          update(nobgb_tolerance_full_model, .~.-(0+Galls|Rhizo)),
          type = "LRT")$`Pr(>Chisq)`[2]
  
  
  nobgb_tolerance_results_final <- results_handling(nobgb_tolerance_results) %>% rowwise() %>%
    mutate(labels2 = case_when(Effect == "Genotype.1" ~ paste("B.", labels2),
                               Effect == "Rhizo.1" ~ paste("C.", labels2),
                               .default = labels2),
           label_border = case_when(Effect %in% c("Genotype.1", "Rhizo.1") ~ Effect,
                                    .default = "Black")) %>%
    as.data.frame()
  
  
  nobgb_tolerance_results_final$labels_x <- c(0.2, 0.25, 0.455, 0.835, NA)
  nobgb_tolerance_results_final$labels_y <- c(0.98,0.97, 0.8875,0.8875, NA)
  
  nobgb_tolerance_VCA <- ggplot(nobgb_tolerance_results_final %>% filter(Effect != "Residuals")) +
    aes(color = Effect, fill = Effect,
        x = mod_VariancePercent,
        pattern = Nema) +
    geom_col_pattern(position = "stack", aes(y = NA),
                     color = "white",
                     linewidth = 1.3,
                     pattern_fill = color_scheme["Nema"],
                     pattern_color = NA) +
    geom_label(data = . %>%
                 #Setting the labels to go in the bar chart vs. beneath
                 filter(Effect %in% c("Genotype", "Genotype:Rhizo")) %>%
                 mutate(labels2 = str_replace(labels2, " \\(", "\\\n\\("),
                        labels2 = str_replace(labels2, "by ", "by\\\n")),
               inherit.aes = FALSE, hjust = 0.5,
               position = "identity",
               aes(label = labels2, y = 1, x = midway)) + 
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x.bottom = element_blank(),
          axis.text.x = element_text(hjust = 0),
          legend.position = "none",
          plot.margin = margin(t=0.2,l=0.4,b=0,r=1, "cm")) +
    xlab("A. Variance explained: host shoot biomass") +
    scale_x_continuous(position = "top", labels = scales::percent) +
    scale_fill_manual(values = color_scheme) +
    scale_pattern_manual(values = c("Infected" = "stripe", "Uninfected" = "none")) +
    theme(plot.margin = margin(t=0.2,l=0.4,b=0,r=0.4, "cm")) +
    coord_cartesian(clip = "off", expand = 0)
  
  nobgb_tolerance_VCA_transition <- ggplot(nobgb_tolerance_results_final %>% filter(Effect != "Residuals") %>% 
                                             #Setting the labels to beneath the bar chart vs. within
                                             filter(!Effect %in% c("Genotype", "Genotype:Rhizo"))
  ) +
    annotate("rect", xmin = 0.025, xmax = 0.49, ymin = 0.05, ymax = 0.8875,
             fill = geno_background, color = color_scheme["Genotype"], linewidth = 2) +
    annotate("rect", xmin = 0.51, xmax = 0.975, ymin = 0.05, ymax = 0.8875,
             fill = rhizo_background, color = color_scheme["Rhizo"], linewidth = 2) +
    geom_segment(inherit.aes = FALSE,
                 aes(x = midway,
                     xend = labels_x,
                     y = 1, 
                     yend = labels_y,
                     color = Effect)) +
    aes(y = labels_y,
        x = labels_x) +
    geom_label(aes(label = labels2, color = label_border),
               label.r = unit(0.5, "lines"),
               hjust = 0.925,
               fill = "white") +
    scale_fill_manual(values = color_scheme, guide = "none") +
    scale_color_manual(values = color_scheme, guide = "none") +
    xlim(c(0,1)) + ylim(c(0,1)) +
    coord_cartesian(expand = 0, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t=0,l=0.4,b=0,r=0.4, "cm"),
          title = element_blank(),
          plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
          panel.background = element_rect(fill = "white", color = "white", linewidth = 2))
  
  # nobgb_tolerance_residuals_model <- update(nobgb_tolerance_full_model,
  # .~.-(Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo)))
  nobgb_tolerance_residuals_model <- update(nobgb_tolerance_full_model, Galls~Below.ground.biomass, family = gaussian)
  
  resids <- residuals(nobgb_tolerance_residuals_model, type = "response")
  length(resids)
  dim(tolerance_data)
  nobgb_tolerance_residual_data <- cbind(tolerance_data, "Residuals" = resids) %>% as.data.frame()
  
  # nobgb_tolerance_vis_model <- glmmTMB(data = nobgb_tolerance_residual_data,
  # formula = Residuals ~ Galls + Genotype + Rhizo + (1|Genotype:Rhizo) + Galls:Genotype + Galls:Rhizo,
  # family = gaussian)
  # nobgb_tolerance_vis_model_GR <- glmmTMB(data = nobgb_tolerance_residual_data %>% mutate(GR = interaction(Genotype, Rhizo)),
  # formula = Residuals ~ Galls + GR + Galls:GR,
  # family = gaussian)
  nobgb_tolerance_vis_model <- update(nobgb_tolerance_full_model,
                                      data = nobgb_tolerance_residual_data,
                                       formula = . ~ Genotype + Rhizo + Residuals + Residuals:Genotype + Residuals:Rhizo)
  
  
  galls_minmax <- plyr::rbind.fill(
    nobgb_tolerance_residual_data %>% group_by(Genotype, Rhizo) %>% summarize(min_galls = min(Galls), max_galls = max(Galls)),
    nobgb_tolerance_residual_data %>% group_by(Rhizo) %>% summarize(min_galls = min(Galls), max_galls = max(Galls)),
    nobgb_tolerance_residual_data %>% group_by(Genotype) %>% summarize(min_galls = min(Galls), max_galls = max(Galls))
  ) %>% mutate(Effect = ifelse(is.na(Rhizo) | is.na(Genotype),"Main" ,"Interaction"))
  
  nobgb_tolerance_viz_rhizo <- emmip(nobgb_tolerance_vis_model,
                                     ~Residuals:Rhizo,
                                     var = "Residuals", 
                                     specs = c("Rhizo", "Residuals"),
                                     type = "response", 
                                     at = list(Residuals = c(seq(-2,4,0.2))),
                                     CIs = TRUE, 
                                     plotit = FALSE, 
                                     rg.limit = 1000000) %>% as.data.frame()
  
  nobgb_tolerance_viz_geno <- emmip(nobgb_tolerance_vis_model,
                                    ~Residuals:Genotype,
                                    var = "Residuals", 
                                    specs = c("Genotype", "Residuals"),
                                    type = "response", 
                                    at = list(Residuals = c(seq(-2,4,0.2))),
                                    CIs = TRUE, 
                                    plotit = FALSE, 
                                    rg.limit = 1000000) %>% as.data.frame()
  
  
  nobgb_tolerance_viz_rhizo_limited <- nobgb_tolerance_viz_rhizo %>%
    left_join(galls_minmax %>% 
                filter(Effect == "Main", !is.na(Rhizo)),
              by = "Rhizo") %>%
    filter(Residuals >= min_galls, Residuals <= max_galls)
  
  nobgb_tolerance_viz_geno_limited <- nobgb_tolerance_viz_geno %>%
    left_join(galls_minmax %>% filter(Effect == "Main", !is.na(Genotype)),
              by = "Genotype") %>%
    filter(Residuals >= min_galls, Residuals <= max_galls)
  
  nobgb_tolerance_max <- ceiling(100 *max(rbind(nobgb_tolerance_viz_rhizo_limited, nobgb_tolerance_viz_geno_limited)$UCL))/100
  nobgb_tolerance_min <- floor(100 *min(rbind(nobgb_tolerance_viz_rhizo_limited, nobgb_tolerance_viz_geno_limited)$LCL))/100
  
  nobgb_tolerance_rhizo <- ggplot(nobgb_tolerance_viz_rhizo_limited) +
    aes(x = as.numeric(Residuals),
        y = yvar,
        ymax = UCL,
        ymin = LCL,
        color = Rhizo, 
        fill = Rhizo) + 
    geom_hline(yintercept = 0, color = "darkgray") + 
    geom_ribbon(alpha = 0.05, color = NA) +
    geom_line() +
    ylim(c(nobgb_tolerance_min, nobgb_tolerance_max)) +
    scale_color_manual(values = rhizo_color_scheme) +
    scale_fill_manual(values = rhizo_color_scheme) +
    ylab("Host Shoot Biomass BLUP") +
    xlab("Gall counts") +
    theme_classic() +
    theme(plot.background = element_rect(fill = NA, color = NA),
          legend.position = "none",
          plot.margin = margin(t = 0.25, b = 1, l = 0.5, r = 0.25, unit = "cm"))
  
  
  nobgb_tolerance_geno <- ggplot(nobgb_tolerance_viz_geno_limited) +
    aes(x = as.numeric(Residuals), y = yvar,
        ymax = UCL, ymin = LCL,
        color = Genotype, fill = Genotype) + 
    geom_hline(yintercept = 0, color = "darkgray") + 
    geom_ribbon(alpha = 0.05, color = NA) +
    geom_line() +
    ylim(c(nobgb_tolerance_min, nobgb_tolerance_max)) +
    scale_color_manual(values = geno_color_scheme) +
    scale_fill_manual(values = geno_color_scheme) +
    ylab("Host Shoot Biomass BLUP") +
    xlab("Gall counts") +
    theme_classic() +
    theme(plot.background = element_rect(fill = NA, color = NA),
          legend.position = "none",
          plot.margin = margin(t = 0.25, b = 1, l = 0.25, r = 0.5, unit = "cm"))
  
  
  nobgb_tolerance_inset_plot <- ggplotGrob(
              plot_grid(nrow = 1, ncol = 4, rel_widths = c(1,12,12,1),
                        blank,
                        nobgb_tolerance_geno,
                        nobgb_tolerance_rhizo,
                        blank)
  )
  
  nobgb_tolerance_model_annotation <- ggplot() + theme_void() +
    theme(plot.background = element_rect(fill = "white", color = "white", linewidth = 2),
          panel.background = element_rect(fill = "white", color = "white", linewidth = 2)) +
    annotation_custom(xmin = 0, xmax = 1, ymin = 0.1, ymax = 0.9,
                      grob = textGrob(x = 0, just = "left",
                                      gp = gpar(fontsize = 8),
                                      paste(
"Variance component analysis model
     Response variable: Shoot biomass
     Fixed effects: Blocking variables, gall counts,
     Random variables: Host genotype, rhizobia strain, host-genotype-by-rhizobia-strain,
                       gall-counts-by-host-genotype, gall-counts-by-rhizoba-strain, gall-counts-by-host-genotype-by-rhizobia-strain\n
     Conditional coefficient of correlation:", round(r.squaredGLMM(nobgb_tolerance_full_model)[2], 4),"
     Marginal coefficient of correlation:", round(r.squaredGLMM(nobgb_tolerance_full_model)[1], 4)
                                      )))
  
  nobgb_tolerance_plot <- plot_grid(ncol = 1, nrow = 4, align = "v", rel_heights = c(1,5,12,4), 
                                    tolerance_title,
                                    nobgb_tolerance_VCA,
                                    nobgb_tolerance_VCA_transition +
                                      annotation_custom(
                                        nobgb_tolerance_inset_plot,
                                        xmin = 0, xmax = 1,
                                        ymin = 0, ymax = 0.875),
                                    nobgb_tolerance_model_annotation
  )
  # nobgb_tolerance_plot
  ggsave(plot = nobgb_tolerance_plot, file = file.path(dir_out, "TolerancePlot_nonbgb.png"),
         width = 11880, height = 8040, unit = "px", dpi = 1200, scale = 1)
}

# nobgb_tolerance_plot
# nobgb_resistance_plot
# nobgb_mutualist

r.squaredGLMM(nobgb_resistance_full_model)
r.squaredGLMM(nobgb_tolerance_full_model)
r.squaredGLMM(nobgb_mutualist_full_model)
