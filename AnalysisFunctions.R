################################################################################
# Author: Addison Dale Buxton-Martin
# 
# Description: This script runs the analysis for ___
# This analysis will import data included in its encasing directory
# and perform analyses as outlined by ___
#
################################################################################

### Defining functions for easier and analyses


check_residuals <- function(model, return_residuals = FALSE, stop_if_nonnormal = FALSE, plot = TRUE, ...){
  #This function will check that the residuals of a model are normally distributed
  #Intended to be used at the beginning of the anlaysis wrapper function
  #Default will stop analysis if there are any important assumptions of an ANOVA test that are not met
  options(contrasts=c("contr.sum", "contr.poly"))
  uniformity_warning <- "! ! ! ! ! ! WARNING: Residuals did not pass the KS test for normality! ! ! ! ! !"
  dispersion_warning <- "! ! ! ! ! ! WARNING: Residuals are not homoscedastic! ! ! ! ! !"
  outlier_warning <- "! ! ! ! ! ! WARNING: Residuals have outliers! ! ! ! ! !"
  stop_warning <- "Model assumptions not met (run with stop_if_nonnormal to proceed)"
  residual_test <- testResiduals(model, plot = FALSE)
  if(plot == TRUE){plot(simulateResiduals(model))}
  if(return_residuals == TRUE){return(residual_test)}
  p_uniformity <- residual_test$uniformity$p.value
  p_dispersion <- residual_test$dispersion$p.value
  p_outlier <- residual_test$outliers$p.value
  if(p_uniformity < 0.05){
    print(uniformity_warning)
    if(stop_if_nonnormal){
      stop(stop_warning)}}
  if(p_dispersion < 0.05){
    print(dispersion_warning)
    if(stop_if_nonnormal){
      stop(stop_warning)}}
  if(p_outlier < 0.05){
    print(outlier_warning)
    if(stop_if_nonnormal){
      stop(stop_warning)}}
}

alias_check <- function(data, formula, print_all_aliases = FALSE, stop_if_alias = FALSE, ...){
  
  alias_warning <- "! ! ! ! ! ! WARNING: Model has aliased coefficients! ! ! ! ! !"
  alias_stop_message <- "Aliased coefficients (run with stop_if_alias as FALSE to proceed)"
  alias_details <- "Showing block aliases only (run with print_all_aliases to see other aliased coefficients; if no table appears an error occurred))"
  
  FE_formula_conversion <- function(string, verbose = TRUE){
    # A function to convert any random effects in a formula into fixed effects
    # Used for alias checking
    out <- paste0(as.character(string))
    #This will remove the encasing "(1|" and ")" from random effects in the formula to make them fixed effects
    out <- gsub("\\(1\\|([a-zA-Z:]+)\\)", "\\1", out)
    #This will remove the encasing "(0+" and ")" from random effects in the formula to make them fixed effects
    out <- gsub("\\(0\\+([a-zA-Z|:()]+)\\)", "\\1", out)
    #This changes the "|" in the above random effects to ":" so they are preserved as interaction effects within a fixed effect framework
    out <- gsub("\\|", ":", out)
    if(verbose){print(paste("Converting random effects to fixed effects: Converted", string, "to", out))}
    return(out)
  }
  #making simple model to check for aliases
  model <- lm(data = data, formula = as.formula(FE_formula_conversion(formula)))
  aliases <- alias(model)
  
  #printing if alises are present
  if(!is.null(aliases$Complete)){
    print(alias_warning)
    print(alias_details)
    try(print(aliases$Complete %>% as.data.frame() %>% select(starts_with("Block")) %>% summarize(across(everything(), ~ sum(. != 0)))))
    try(print(aliases$Complete %>% as.data.frame() %>% t() %>% as.data.frame() %>% select(starts_with("Block")) %>% summarize(across(everything(), ~ sum(. != 0)))))
    if(print_all_aliases){print(aliases)}
    if(stop_if_alias == TRUE){stop(alias_stop_message)}
  }
}


FE_analysis <- function(formula,
                        df = my_data,
                        family_distribution = nbinom1,
                        return_model = FALSE,
                        singular.ok = FALSE,
                        block = block_formula,
                        ...){
  #Setting proper statistical options
  options(contrasts=c("contr.sum", "contr.poly"))
  
  #Argument handling: ensuring that only the arguments in FE_analysis, check_residuals,
  #and alias_check are being passed to the FE_analysis function
  argg <- names(c(as.list(environment()), list(...)))
  total_formals <- c(names(formals(FE_analysis)), names(formals(check_residuals)), names(formals(alias_check)))
  if(sum(argg %in% total_formals) != length(argg)){
    print(argg[which(!argg %in% total_formals)])
    stop("Unusable argument used")}
  
  #Dealing with formula argument and the block part of the argument
  if(block == TRUE){
    form <- paste(formula, block_formula)}
  else if(block == FALSE){
    form <- paste(formula)}
  else {
    form <- paste(formula, block)}
  
  #Fitting data to model
  model <- glmmTMB(data = df, family = family_distribution, formula = as.formula(form))
  
  #Checking model fit and for aliased variables
  print("Checking model fit")
  check_residuals(model = model, ...)
  alias_check(data = df, formula = form, ...)
  
  #Printing anova and returning model if required
  anova_results <- Anova(type = 3, model, singular.ok = singular.ok)
  anova_results <- anova_results %>% mutate(sig = case_when(`Pr(>Chisq)` < 0.0001 ~ "****",
                                                            `Pr(>Chisq)` < 0.001 ~ "***",
                                                            `Pr(>Chisq)` < 0.01 ~ "**",
                                                            `Pr(>Chisq)` < 0.05 ~ "*",
                                                            `Pr(>Chisq)` < 0.1 ~ ".",
                                                            .default = ""))
  print(anova_results)
  if(return_model){return(model)}else(return(anova_results %>% as.data.frame()))
}

variance_to_table <- function(varcor){
  return(varcor$cond %>% 
           as_tibble() %>% t() %>% as.data.frame() %>% 
           rename("Variance" = "V1") %>%
           rownames_to_column("Effect"))}

generate_formula_table <- function(trait, volume = FALSE, galls = FALSE, mod_trait = FALSE){
  if(!(trait %in% colnames(my_data))){stop("Error in generate_formula_table: trait must be a column name in the data table (my_data)")}
  if(!(volume == TRUE | volume == FALSE | volume == "scale")){stop("Error in generate_formula_table: volume must be TRUE, FALSE, or 'scale'")}
  if(!(galls == TRUE | galls == FALSE | galls == "scale")){stop("Error in generate_formula_table: galls must be TRUE, FALSE, or 'scale'")}
  if(!(mod_trait == FALSE | is.character(mod_trait))){stop("Error in generate_formula_table: mod_trait must be FALSE or a string")}
  
  nema_table <- rbind(
    c("Effect" = "Nema:Genotype:Rhizo",
      "Full" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo) + (1|Nema:Genotype:Rhizo)",
      "Reduced" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)"),
    c("Effect" = "Nema:Rhizo", 
      "Full" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)",
      "Reduced" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype)"),
    c("Effect" = "Nema:Genotype",
      "Full" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)",
      "Reduced" = "Trait ~ Nema +  (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Rhizo)"),
    c("Effect" = "Genotype:Rhizo",
      "Full" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)",
      "Reduced" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)"),
    c("Effect" = "Genotype",
      "Full" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo)",
      "Reduced" = "Trait ~ Nema + (1|Rhizo)"),
    c("Effect" = "Rhizo",
      "Full" = "Trait ~ Nema + (1|Genotype) + (1|Rhizo)",
      "Reduced" = "Trait ~ Nema + (1|Genotype)"))
  gall_table <- rbind(
    c("Effect" = "Genotype.Rhizo.1", #Genotype:Rhizo:Galls
      "Full" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo)",
      "Reduced" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)"),
    c("Effect" = "Rhizo.1", #Rhizo:Galls
      "Full" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
      "Reduced" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype)"),
    c("Effect" = "Genotype.1", #Genotype:Galls
      "Full" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
      "Reduced" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Rhizo)"),
    c("Effect" = "Genotype.Rhizo", #Genotype:Galls
      "Full" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
      "Reduced" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)"),
    c("Effect" = "Genotype",
      "Full" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo)",
      "Reduced" = "Trait ~ Galls + (1|Rhizo)"),
    c("Effect" = "Rhizo",
      "Full" = "Trait ~ Galls + (1|Genotype) + (1|Rhizo)",
      "Reduced" = "Trait ~ Galls + (1|Genotype)"))
  nema_vol_table <- rbind(
    c("Effect" = "Nema:Genotype:Rhizo",
      "Full" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo) + (1|Nema:Genotype:Rhizo)",
      "Reduced" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)"),
    c("Effect" = "Nema:Rhizo", 
      "Full" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)",
      "Reduced" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype)"),
    c("Effect" = "Nema:Genotype",
      "Full" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)",
      "Reduced" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Rhizo)"),
    c("Effect" = "Genotype:Rhizo",
      "Full" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)",
      "Reduced" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Nema:Genotype) + (1|Nema:Rhizo)"),
    c("Effect" = "Genotype",
      "Full" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo)",
      "Reduced" = "Trait ~ Nema + Volume.mm3 + (1|Rhizo)"),
    c("Effect" = "Rhizo",
      "Full" = "Trait ~ Nema + Volume.mm3 + (1|Genotype) + (1|Rhizo)",
      "Reduced" = "Trait ~ Nema + Volume.mm3 + (1|Genotype)"))
  gall_vol_table <- rbind(
    c("Effect" = "Genotype.Rhizo.1", #Genotype:Rhizo:Galls
      "Full" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo) + (0+Galls|Genotype:Rhizo)",
      "Reduced" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)"),
    c("Effect" = "Rhizo.1", #Rhizo:Galls
      "Full" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
      "Reduced" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype)"),
    c("Effect" = "Genotype.1", #Genotype:Galls
      "Full" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
      "Reduced" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Rhizo)"),
    c("Effect" = "Genotype.Rhizo", #Genotype:Galls
      "Full" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (1|Genotype:Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)",
      "Reduced" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo) + (0+Galls|Genotype) + (0+Galls|Rhizo)"),
    c("Effect" = "Genotype",
      "Full" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo)",
      "Reduced" = "Trait ~ Galls + Volume.mm3 + (1|Rhizo)"),
    c("Effect" = "Rhizo",
      "Full" = "Trait ~ Galls + Volume.mm3 + (1|Genotype) + (1|Rhizo)",
      "Reduced" = "Trait ~ Galls + Volume.mm3 + (1|Genotype)"))
  
  if(galls == FALSE){
    if(volume == FALSE){
      table <- nema_table} 
    else {
        table <- nema_vol_table}}
  else {
    if(volume == FALSE){
      table <- gall_table}
    else {
        table <- gall_vol_table}}
  
  if(volume == "scale"){table <- gsub("Volume.mm3", "scale(Volume.mm3)", table)}
  if(galls == "scale"){table <- gsub("Galls", "scale(Galls)", table)}
  if(is.character(mod_trait)){table <- gsub("Trait", paste0(mod_trait,"(Trait)"), table)}
  table <- gsub("Trait", trait, table)
  return(table)
  }

RE_analysis <- function(formula_table,
                        fullest_formula,
                        df = my_data,
                        family_distribution = nbinom1,
                        return_model = FALSE,
                        singular.ok = FALSE,
                        block = block_formula,
                        recursive = TRUE, ...){
  ### This formula takes as its input a formula_table, where the first column
  # is called "Effect" and has all of the Effect that you are interested in.
  # The effects must correspond to all of the variables in the fullest_formula
  # variable which is a string representing the fullest formula you are using
  # in the analysis. The second column of the formula table is called "Full"
  # and contains the full formula used in a likelihood ratio test to test whether
  # a variance component is different from zero. The third column is called "Reduced"
  # and is the reduced formula for the likelihood ratio test.
  
  #Setting proper statistical options
  options(contrasts=c("contr.sum", "contr.poly"))
  
  #Argument handling: ensuring that only the arguments in FE_analysis, check_residuals,
  #and alias_check are being passed to the FE_analysis function
  argg <- names(c(as.list(environment()), list(...)))
  total_formals <- c(names(formals(RE_analysis)), names(formals(check_residuals)), names(formals(alias_check)))
  if(sum(argg %in% total_formals) != length(argg)){
    print(argg[which(!argg %in% total_formals)])
    stop("Unusable argument used")}
  
  #Dealing with formula argument and the block part of the argument
  formula_table <- formula_table %>% as.data.frame() %>% rowwise() %>%
    mutate(Full = if_else(block == TRUE, paste(Full, block_formula), if_else(block == FALSE, Full, paste(Full, block))),
           Reduced = if_else(block == TRUE, paste(Reduced, block_formula), if_else(block == FALSE, Reduced, paste(Reduced, block))))
  
  if(block == TRUE){
    form <- paste(fullest_formula, block_formula)}
  else if(block == FALSE){
    form <- paste(fullest_formula)}
  else {
    form <- paste(fullest_formula, block)}
  
  # Model for variances
  fullest_model <- glmmTMB(data = df, family = family_distribution, formula = as.formula(form))
  Fixed_effect_table <- Anova(type = 3, fullest_model) %>% rownames_to_column(var = "Effect") %>% 
    rename("Chi Df" = "Df", "pvalue" = "Pr(>Chisq)") %>% as.data.frame()
  print(colnames(Fixed_effect_table))
  
  #Checking model fit and for aliased variables
  print("Checking model fit")
  check_residuals(model = fullest_model, ...)
  alias_check(data = df, formula = form, ...)
  
  if(return_model){return(fullest_model)}else{
  #Variances
  print("Getting variances")
  varcor_results <- variance_to_table(VarCorr(fullest_model))
  
  #Defining function for getting p values of likelihood estimates
  analysis <- function(full_formula, reduced_formula, df, family_distribution, singular.ok){
    
    #Fitting data to model
    full_model <- glmmTMB(data = df, family = family_distribution, formula = as.formula(full_formula))
    red_model <- glmmTMB(data = df, family = family_distribution, formula = as.formula(reduced_formula))
    results <- anova(full_model, red_model, test = "LRT") %>% as.data.frame()
    returned_results <- c("Chisq" = results[2, "Chisq"], "Chidf" = results[2, "Chi Df"], "pvalue" = results[2,"Pr(>Chisq)"])
    returned_results <- returned_results #%>% as.data.frame()
    return(returned_results)
  }
  
    #Updating table
  lrt_results <- apply(formula_table, 1, function(row) {
    res <- analysis(df = df, family_distribution = family_distribution, singular.ok = singular.ok,
                    full_formula = row[2], reduced_formula = row[3])
    return(c(row[[1]], res[[1]], res[[2]], res[[3]]))
  })
  lrt_results <- t(lrt_results) %>% as.data.frame()
  names(lrt_results) <- c("Effect", "Chisq", "Chi Df", "pvalue")
  updated_table <- full_join(as.data.frame(formula_table), as.data.frame(lrt_results), by = "Effect")
  final_results <- full_join(updated_table, varcor_results, by = "Effect")
  print(Fixed_effect_table)
  final_results <- plyr::rbind.fill(final_results %>% as.data.frame(), Fixed_effect_table %>% as.data.frame())
  final_results <- final_results %>% mutate(sig = case_when(as.numeric(pvalue) < 0.0001 ~ "****",
                                                            as.numeric(pvalue) < 0.001 ~ "***",
                                                            as.numeric(pvalue) < 0.01 ~ "**",
                                                            as.numeric(pvalue) < 0.05 ~ "*",
                                                            as.numeric(pvalue) < 0.1 ~ ".",
                                                            .default = ""))
  return(final_results[,c(1,7,4,5,6,8,2,3)])
  }
}