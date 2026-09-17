#################################################################################################
# General information 
#################################################################################################

# Author: Marie de Bakker
# Date: 11/09/2024
# Project: iPREVENT-CVD

# Information: functions

library(ggplot2)

###################
#### function 1 ###
###################

predictions_for_plot <- function(model, dataset, variable_yaxis, variable_xaxis = "age"){
  
  # preparations  ----------------------------------------------------------------------------------------------------------
  # create unified columns to reference y-axis and x-axis variables
  dataset$variable_yaxis <- dataset[, variable_yaxis]
  dataset$variable_xaxis <- dataset[, variable_xaxis]
  
  # extract fixed effects part of the model formula
  fixed_effects <- paste0(formula(model))[3]
  model_formula_fixed <- paste0("~ ", fixed_effects)
  
  # define x-axis range for plotting (excluding extreme outliers)
  xaxis_lower <- quantile(dataset[!is.na(dataset$variable_yaxis), "variable_xaxis"], 0.005, na.rm = T) %>% as.numeric() %>% round(digits = 0)
  xaxis_upper <- quantile(dataset[!is.na(dataset$variable_yaxis), "variable_xaxis"], 0.995, na.rm = T) %>% as.numeric() %>% round(digits = 0)
  
  # create df_new  ----------------------------------------------------------------------------------------------------------
  df_new <- with(dataset, expand.grid(variable_xaxis = seq(xaxis_lower, xaxis_upper, length.out = 250),
                                      sex = levels(sex),
                                      med_antihypertensive = levels(med_antihypertensive),
                                      med_lipidlowering = levels(med_lipidlowering),
                                      med_antidiabetic = levels(med_antidiabetic),
                                      cohort_name = levels(cohort_name)))
  
  # set contrasts for cohort variable
  contrasts(df_new$cohort_name) <- contr.sum(levels(dataset$cohort_name))
  
  # rename variable_xaxis to its original name to make model.matrix work
  colnames(df_new)[1] <- variable_xaxis
  
  # calculate predictions ----------------------------------------------------------------------------------------------------------
  
  # for numeric outcomes ------------------------
  if(is.numeric(dataset$variable_yaxis)){
    
    # build model matrix
    X <- model.matrix(as.formula(model_formula_fixed), data = df_new)
    
    # set dummy cohort variables to 0 for global effect
    X[, grepl('^cohort', colnames(X))] <- 0
    
    # get fixed effect coefficients and covariance matrix
    betas <- fixef(model)
    V <- vcov(model)
    
    # compute predictions on log2 scale
    df_new$log2_pred <- c(X %*% betas) 
    df_new$ses <- sqrt(diag(X %*% V %*% t(X)))
    df_new$log2_lowerlimit <- df_new$log2_pred - 1.96*df_new$ses
    df_new$log2_upperlimit <- df_new$log2_pred + 1.96*df_new$ses
    
    # back-transform to original scale unless already standardized and log2-tranformed
    if (!grepl("^zlog2", variable_yaxis)) {
      df_new$pred <- 2^(df_new$log2_pred)
      df_new$lowerlimit <- 2^(df_new$log2_lowerlimit)
      df_new$upperlimit <- 2^(df_new$log2_upperlimit)
    } else if (grepl("^zlog2", variable_yaxis)){
      df_new$pred <- df_new$log2_pred
      df_new$lowerlimit <- df_new$log2_lowerlimit
      df_new$upperlimit <- df_new$log2_upperlimit
    }
    
    # for factor outcomes ------------------------
  } else if(is.factor(dataset$variable_yaxis)){
    
    # create df_new  ----------------------------------------------------------------------------------------------------------
    df_new <- with(dataset, expand.grid(variable_xaxis = seq(xaxis_lower, xaxis_upper, length.out = 100),
                                        sex = levels(sex),
                                        med_antihypertensive = levels(med_antihypertensive),
                                        med_lipidlowering = levels(med_lipidlowering),
                                        med_antidiabetic = levels(med_antidiabetic),
                                        cohort_name1 = 0,
                                        cohort_name2 = 0,
                                        cohort_name3 = 0))
    
    # rename variable_xaxis back to original name
    colnames(df_new)[1] <- variable_xaxis
    
    # build model matrix
    termsX <- delete.response(model$Terms$termsX)
    mfX <- model.frame(termsX, df_new, 
                       xlev = .getXlevels(termsX, model$model_frames$mfX))
    X <- model.matrix(termsX, mfX)
    
    # obtain marginal coefficients
    mcoefs <- marginal_coefs(model, std_errors = TRUE, seed = 1, cores = 13)
    betas <- mcoefs$betas
    var_betas <- mcoefs$var_betas
    
    # get predictions on logit scale and transform to probability scale
    df_new$exp_pred <- c(X %*% betas)
    df_new$ses <- sqrt(diag(X %*% var_betas %*% t(X)))
    df_new$exp_lowerlimit <- df_new$exp_pred + (qnorm((1 - 0.95) / 2) * df_new$ses)
    df_new$exp_upperlimit <- df_new$exp_pred + (qnorm((1 + 0.95) / 2) * df_new$ses)
    
    # transform logits to probabilities (percent)
    df_new$pred <- (exp(df_new$exp_pred) / (1 + exp(df_new$exp_pred)))*100
    df_new$lowerlimit <- (exp(df_new$exp_lowerlimit) / (1 + exp(df_new$exp_lowerlimit)))*100
    df_new$upperlimit <- (exp(df_new$exp_upperlimit) / (1 + exp(df_new$exp_upperlimit)))*100
  }
  
  # return dataframe with predictions  -------------------------------------------------------------------------------------------------
  return(df_new)
}

###################
#### function 2 ###
###################

predictions_for_plot_stratified <- function(model, dataset, variable_yaxis, variable_xaxis = "age"){
  
  # preparations  ----------------------------------------------------------------------------------------------------------
  # create unified columns to reference y-axis and x-axis variables
  dataset$variable_yaxis <- dataset[, variable_yaxis]
  dataset$variable_xaxis <- dataset[, variable_xaxis]
  
  # extract fixed effects part of the model formula
  fixed_effects <- paste0(formula(model))[3]
  model_formula_fixed <- paste0("~ ", fixed_effects)
  
  # define x-axis range for plotting (excluding extreme outliers)
  xaxis_lower <- quantile(dataset[!is.na(dataset$variable_yaxis), "variable_xaxis"], 0.005, na.rm = T) %>% as.numeric() %>% round(digits = 0)
  xaxis_upper <- quantile(dataset[!is.na(dataset$variable_yaxis), "variable_xaxis"], 0.995, na.rm = T) %>% as.numeric() %>% round(digits = 0)
  
  # create df_new  ----------------------------------------------------------------------------------------------------------
  df_new <- with(dataset, expand.grid(variable_xaxis = seq(xaxis_lower, xaxis_upper, length.out = 250),
                                      sex = levels(sex),
                                      med_antihypertensive = levels(med_antihypertensive),
                                      med_lipidlowering = levels(med_lipidlowering),
                                      med_antidiabetic = levels(med_antidiabetic)))
  
  # rename variable_xaxis to its original name to make model.matrix work
  colnames(df_new)[1] <- variable_xaxis
  
  # calculate predictions ----------------------------------------------------------------------------------------------------------
  
  # for numeric outcomes ------------------------
  if(is.numeric(dataset$variable_yaxis)){
    # build model matrix
    X <- model.matrix(as.formula(model_formula_fixed), data = df_new)
    
    # get fixed effect coefficients and covariance matrix
    betas <- fixef(model)
    V <- vcov(model)
    
    # compute predictions on log2 scale
    df_new$log2_pred <- c(X %*% betas) 
    df_new$ses <- sqrt(diag(X %*% V %*% t(X)))
    df_new$log2_lowerlimit <- df_new$log2_pred - 1.96*df_new$ses
    df_new$log2_upperlimit <- df_new$log2_pred + 1.96*df_new$ses
    
    # back-transform to original scale unless already standardized and log2-tranformed
    if (!grepl("^zlog2", variable_yaxis)) {
      df_new$pred <- 2^(df_new$log2_pred)
      df_new$lowerlimit <- 2^(df_new$log2_lowerlimit)
      df_new$upperlimit <- 2^(df_new$log2_upperlimit)
    } else if (grepl("^zlog2", variable_yaxis)){
      df_new$pred <- df_new$log2_pred
      df_new$lowerlimit <- df_new$log2_lowerlimit
      df_new$upperlimit <- df_new$log2_upperlimit
    }
    
    # for factor outcomes ------------------------
  } else if(is.factor(dataset$variable_yaxis)){
    
    # get predictions on logit scale and transform to probability scale
    df_new <- effectPlotData(model, df_new, marginal = T, cores = 12) %>% rename(exp_pred = pred,
                                                                     exp_lowerlimit = low,
                                                                     exp_upperlimit = upp)
    
    # transform logits to probabilities (percent)
    df_new$pred <- (exp(df_new$exp_pred) / (1 + exp(df_new$exp_pred)))*100
    df_new$lowerlimit <- (exp(df_new$exp_lowerlimit) / (1 + exp(df_new$exp_lowerlimit)))*100
    df_new$upperlimit <- (exp(df_new$exp_upperlimit) / (1 + exp(df_new$exp_upperlimit)))*100
  }
  
  # return dataframe with predictions  -------------------------------------------------------------------------------------------------
  return(df_new)
}


###################
#### function 3 ###
###################

trajectory_plot_age <- function(df_new,
                                plot_antihypertensive = "No", plot_lipidlowering = "No", plot_antidiabetic = "No",
                                colors = c("#DD3C51", "#1F6683"),
                                limits_xaxis = c(10, 100), limits_yaxis = NULL,
                                plot_title = "", label_xaxis = "Age (years)", label_yaxis = "",
                                legend_labels = c("Women", "Men"),
                                y.ticks_multiplier = 2.5, y.at_multiplier = 5, legend_positions = c(1.01, 1.0)){
  
  # determine limits of xaxis
  if(is.null(limits_xaxis)){
    xlim_lower <- round(min(df_new$age, na.rm = T),1)
    xlim_upper <- round(max(df_new$age, na.rm = T),1)
  } else {
    xlim_lower <- limits_xaxis[1]
    xlim_upper <- limits_xaxis[2]
  }
  
  # select 'subjects' that should be plotted
  df_plot <- subset(df_new, med_antihypertensive == plot_antihypertensive & med_lipidlowering == plot_lipidlowering &  med_antidiabetic == plot_antidiabetic &
                      age >= xlim_lower & age <= xlim_upper)
  
  # preparations for plot  ----------------------------------------------------------------
  
  # set title and labels for x-axis/y-axis
  label_yaxis = label_yaxis
  label_xaxis = label_xaxis
  plot_title = plot_title
  legend_labels = legend_labels
  
  # determine limits of yaxis
  if(is.null(limits_yaxis)){
    ylim_upper <- ceiling(max(df_plot$pred))
    ylim_lower <- floor(min(df_plot$pred))
  } else {
    ylim_lower <- limits_yaxis[1]
    ylim_upper <- limits_yaxis[2]
  }
  
  # set ticks for xaxis
  x.at = seq(xlim_lower, xlim_upper, 10) # text
  x.ticks = seq(xlim_lower, xlim_upper, 5) # ticks
  
  # set ticks for yaxis 
  y.at = seq(ylim_lower,ylim_upper, y.at_multiplier) # text
  y.ticks = seq(ylim_lower, ylim_upper, y.ticks_multiplier) # ticks
  
  # plit dataset for sex
  df_women <- subset(df_plot, sex == "Female")
  df_men <- subset(df_plot, sex == "Male")
  
  # plot  ----------------------------------------------------------------
  
  # plot (ggplot)
  plot <- ggplot() +
    # event_1
    geom_ribbon(data = df_women, aes(x=age, ymin=lowerlimit, ymax=upperlimit), alpha=0.1, fill=colors[1]) +
    geom_line(data=df_women, aes(x=age, y=pred), color=colors[1], linewidth = 1.25) +
    annotate("text", 
             label=legend_labels[1], 
             x=xlim_lower+1, 
             y=ylim_upper*legend_positions[1],
             col=colors[1], 
             hjust="left", 
             size = 5) +
    # event_2
    geom_ribbon(data = df_men, aes(x=age, ymin=lowerlimit, ymax=upperlimit), alpha=0.1, fill=colors[2]) +
    geom_line(data=df_men, aes(x=age, y=pred), color=colors[2], linewidth = 1.25) +
    annotate("text", 
             label=legend_labels[2], 
             x=xlim_lower+1, 
             y=ylim_upper*legend_positions[2],
             col=colors[2], 
             hjust="left", 
             size = 5) +
    # titles and grid
    labs(x=label_xaxis, y=make_label(label_yaxis)) +
    ggtitle(plot_title) +
    theme(panel.background = element_rect(fill = "white", linewidth =  2, linetype = "solid"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey"), 
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size=16),
          axis.text.y = element_text(size=16),
          axis.text.x = element_text(size=16),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    scale_x_continuous(minor_breaks = x.ticks, breaks = c(x.at), expand=expansion(mult=c(0.05,0.05)))  +
    scale_y_continuous(minor_breaks = y.ticks, breaks = c(y.at), expand=expansion(mult=c(0.15,0.15))) +
    coord_cartesian(ylim=c(ylim_lower, ylim_upper), xlim=c(xlim_lower, xlim_upper))
  
  return(plot)
  
}

###################
#### function 4 ###
###################

plot_residuals <- function(sample_size = 0.3, fitted_model, level = 0){
  # Define the fraction of points to sample (e.g., 30% of data)
  set.seed(123)  # For reproducibility
  sample_size <- sample_size  # Adjust this fraction as needed
  model <- fitted_model
  
  # model class
  model_class <- class(model)
  
  # Extract fitted values and residuals based on model class
  if (inherits(model, "lme")) {  # nlme package
    fitted_vals <- fitted(model, level = level)
    resid_vals <- resid(model, type = "p", level = level)
    model_data <- getData(model)  # Retrieve model's actual data
  } else if (inherits(model, "merMod")) {  # lme4 package
    fitted_vals <- fitted(model)
    resid_vals <- residuals(model, type = "pearson")
    model_data <- model@frame  # Extract model data
  } else if (inherits(model, "MixMod")) {  # GLMMadaptive package
    fitted_vals <- fitted(model, level = level)  # Fitted values
    if(level == 0) { resid_vals <- resid(model, type = "mean_subject")
    } else if (level == 1){resid_vals <- resid(model, type = "subject_specific")}
    model_data <- model$data  # Retrieve model data
  } else {
    stop("Unsupported model type: ", paste(model_class, collapse = ", "), 
         ". Provide a model from 'nlme' (lme), 'lme4' (merMod), or 'GLMMadaptive' (MixMod).")
  }
  
  # Ensure they have the same length (handle possible missing values)
  valid_indices <- complete.cases(fitted_vals, resid_vals)  # Find complete cases
  fitted_vals <- fitted_vals[valid_indices]
  resid_vals <- resid_vals[valid_indices]
  
  # Extract grouping variable (sex)
  if (inherits(model, "lme")) {  # nlme package
    sex_vals <- model_data$sex[valid_indices]  # Align sex variable correctly
  } else if (inherits(model, "MixMod")) {  # GLMMadaptive package
    sex_vals <- model_data$sex[names(resid_vals) %>% as.numeric]
  } 
  
  # Create a data frame with valid values
  resid_data <- data.frame(fitted = fitted_vals, resid = resid_vals, sex = sex_vals)
  
  # Sample a subset of data
  sample_indices <- sample(nrow(resid_data), size = floor(sample_size * nrow(resid_data)))
  resid_sample <- resid_data[sample_indices, ]
  
  # Split data into males and females
  resid_male <- subset(resid_sample, sex == "Male")
  resid_female <- subset(resid_sample, sex == "Female")
  
  # Determine common axis limits
  x_limits <- range(resid_sample$fitted)  # Common x-axis (fitted values)
  y_limits <- range(resid_sample$resid)   # Common y-axis (residuals)
  
  # Set up side-by-side plotting
  par(mfrow = c(1, 2))  # 1 row, 2 columns
  
  # Plot for Males
  plot(resid_male$fitted, resid_male$resid, col = "#1F6683", pch = 16,
       main = "Residuals vs. Fitted (Males)",
       xlab = "Fitted Values", ylab = "Standardized Residuals",
       xlim = x_limits, ylim = y_limits)
  lines(lowess(resid_male$fitted, resid_male$resid), col = "darkgrey", lwd = 3)
  
  # Plot for Females
  plot(resid_female$fitted, resid_female$resid, col = "#DD3C51", pch = 16,
       main = "Residuals vs. Fitted (Females)",
       xlab = "Fitted Values", ylab = "Standardized Residuals",
       xlim = x_limits, ylim = y_limits)
  lines(lowess(resid_female$fitted, resid_female$resid), col = "darkgrey", lwd = 3)
  
  # Reset plotting layout
  par(mfrow = c(1, 1))
  
}


###################
#### function 5 ###
###################

age_women_exceeding_men <- function(df_new, limits_xaxis = NULL,
                                    plot_antihypertensive = "No", 
                                    plot_lipidlowering = "No",
                                    plot_antidiabetic = "No",
                                    variable_yaxis = NULL){
  
  # set seed for reproducibility in random sampling
  set.seed(123)
  number_of_samples <- 10000  # Number of random samples for bootstrapping
  
  # define a helper function that generates samples from a normal distribution
  # with given mean (log2_pred) and standard deviation (ses)
  sampler <- function(x, y) rnorm(number_of_samples, mean = x, sd = y)
  
  # define x-axis limits either from input or inferred from data
  # (to make sure we select the correct data for this function to match original plot)
  if(is.null(limits_xaxis)){
    xlim_lower <- round(min(df_new$age, na.rm = TRUE), 1)
    xlim_upper <- round(max(df_new$age, na.rm = TRUE), 1)
  } else {
    xlim_lower <- limits_xaxis[1]
    xlim_upper <- limits_xaxis[2]
  }
  
  # filter dataset for selected medication status and age range
  # (to make sure we select the correct data for this function to match original plot)
  # since predictions are the same for each cohort ('average' cohort effect), select data
  # for only one cohort.
  if(is.null(variable_yaxis)){
    df_plot <- subset(df_new, 
                      med_antihypertensive == plot_antihypertensive &
                        med_lipidlowering == plot_lipidlowering &
                        med_antidiabetic == plot_antidiabetic &
                        age >= xlim_lower & age <= xlim_upper & 
                        cohort_name == levels(cohort_name)[1])
  } else if(variable_yaxis == "current_smoker"){
    df_plot <- subset(df_new, 
                      med_antihypertensive == plot_antihypertensive &
                        med_lipidlowering == plot_lipidlowering &
                        med_antidiabetic == plot_antidiabetic &
                        age >= xlim_lower & age <= xlim_upper)
  }

  
  df_plot_f <- subset(df_plot, sex == "Female")
  df_plot_m <- subset(df_plot, sex == "Male")
  
  # Sort by age to align indices
  df_plot_f <- df_plot_f[order(df_plot_f$age), ]
  df_plot_m <- df_plot_m[order(df_plot_m$age), ]
  ages <- df_plot_f$age  # ordered list of ages
  
  # generate random samples for each age based on log2 prediction and SE  (matrix: age x sample)
  # this approximates the uncertainty in prediction
  samples_f <- mapply(sampler, df_plot_f$log2_pred, df_plot_f$ses)
  samples_m <- mapply(sampler, df_plot_m$log2_pred, df_plot_m$ses)
  
  # transpose so we get samples in rows and ages in columns
  samples_f <- t(samples_f)  # [age x sample]
  samples_m <- t(samples_m)
  
  # now loop through each sample and find the first age where females > males
  ages_exceed <- numeric(number_of_samples)
  
  for(i in 1:number_of_samples){
    diff <- samples_f[, i] - samples_m[, i]
    exceed_idx <- which(diff > 0)
    
    if(length(exceed_idx) > 0){
      ages_exceed[i] <- ages[min(exceed_idx)]
    } else {
      ages_exceed[i] <- NA  # in case females never exceed males
    }
  }
  
  # remove NAs if any
  ages_exceed <- na.omit(ages_exceed)
  
  # Get CI
  median_age <- median(ages_exceed)
  lowerlimit_age <- quantile(ages_exceed, 0.025)
  upperlimit_age <- quantile(ages_exceed, 0.975)
  
  
  # print result
  print(paste0("Females exceed males at median age: ",
               round(median_age, 1),
               " (95% CI: ",
               round(lowerlimit_age, 1),
               " to ",
               round(upperlimit_age, 1), ")"))
}


###################
#### function 6 ###
###################

make_label <- function(label) {
  
  # replace spaces with plotmath spacing
  label <- gsub(" ", "~", label)
  
  # subscripts for logX
  label <- gsub("log([0-9]+)", "log[\\1]", label)
  
  # superscripts for m2, m3, etc.
  label <- gsub("m([0-9]+)", "m^\\1", label)
  
  # return as expression
  as.expression(parse(text = label))
}



###################
#### function 5 ###
###################

age_men_exceeding_women <- function(df_new, limits_xaxis = NULL,
                                    plot_antihypertensive = "No", 
                                    plot_lipidlowering = "No",
                                    plot_antidiabetic = "No",
                                    variable_yaxis = NULL){
  
  # set seed for reproducibility in random sampling
  set.seed(123)
  number_of_samples <- 10000  # Number of random samples for bootstrapping
  
  # define a helper function that generates samples from a normal distribution
  # with given mean (log2_pred) and standard deviation (ses)
  sampler <- function(x, y) rnorm(number_of_samples, mean = x, sd = y)
  
  # define x-axis limits either from input or inferred from data
  # (to make sure we select the correct data for this function to match original plot)
  if(is.null(limits_xaxis)){
    xlim_lower <- round(min(df_new$age, na.rm = TRUE), 1)
    xlim_upper <- round(max(df_new$age, na.rm = TRUE), 1)
  } else {
    xlim_lower <- limits_xaxis[1]
    xlim_upper <- limits_xaxis[2]
  }
  
  # filter dataset for selected medication status and age range
  # (to make sure we select the correct data for this function to match original plot)
  # since predictions are the same for each cohort ('average' cohort effect), select data
  # for only one cohort.
  if(is.null(variable_yaxis)){
    df_plot <- subset(df_new, 
                      med_antihypertensive == plot_antihypertensive &
                        med_lipidlowering == plot_lipidlowering &
                        med_antidiabetic == plot_antidiabetic &
                        age >= xlim_lower & age <= xlim_upper & 
                        cohort_name == levels(cohort_name)[1])
  } else if(variable_yaxis == "current_smoker"){
    df_plot <- subset(df_new, 
                      med_antihypertensive == plot_antihypertensive &
                        med_lipidlowering == plot_lipidlowering &
                        med_antidiabetic == plot_antidiabetic &
                        age >= xlim_lower & age <= xlim_upper)
  }
  
  
  df_plot_f <- subset(df_plot, sex == "Female")
  df_plot_m <- subset(df_plot, sex == "Male")
  
  # Sort by age to align indices
  df_plot_f <- df_plot_f[order(df_plot_f$age), ]
  df_plot_m <- df_plot_m[order(df_plot_m$age), ]
  ages <- df_plot_f$age  # ordered list of ages
  
  # generate random samples for each age based on log2 prediction and SE  (matrix: age x sample)
  # this approximates the uncertainty in prediction
  samples_f <- mapply(sampler, df_plot_f$log2_pred, df_plot_f$ses)
  samples_m <- mapply(sampler, df_plot_m$log2_pred, df_plot_m$ses)
  
  # transpose so we get samples in rows and ages in columns
  samples_f <- t(samples_f)  # [age x sample]
  samples_m <- t(samples_m)
  
  # now loop through each sample and find the first age where females > males
  ages_exceed <- numeric(number_of_samples)
  
  for(i in 1:number_of_samples){
    diff <- samples_m[, i] - samples_f[, i]
    exceed_idx <- which(diff > 0)
    
    if(length(exceed_idx) > 0){
      ages_exceed[i] <- ages[min(exceed_idx)]
    } else {
      ages_exceed[i] <- NA  # in case females never exceed males
    }
  }
  
  # remove NAs if any
  ages_exceed <- na.omit(ages_exceed)
  
  # Get CI
  median_age <- median(ages_exceed)
  lowerlimit_age <- quantile(ages_exceed, 0.025)
  upperlimit_age <- quantile(ages_exceed, 0.975)
  
  
  # print result
  print(paste0("Males exceed females at median age: ",
               round(median_age, 1),
               " (95% CI: ",
               round(lowerlimit_age, 1),
               " to ",
               round(upperlimit_age, 1), ")"))
}



