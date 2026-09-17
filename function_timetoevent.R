#################################################################################################
# General information 
#################################################################################################

# Author: Marie de Bakker
# Date: 11/09/2024
# Project: iPREVENT-CVD

# Information: functions to predict effects from mixed effect models and plot trajectories

###################
#### function 1 ###
###################

predictions_for_plot <- function(model, dataset, variable_yaxis, variable_xaxis = "time_to_event", plot_age = NULL){
  
  # preparations  ----------------------------------------------------------------------------------------------------------
  # create unified columns to reference y-axis and x-axis variables
  dataset$variable_yaxis <- dataset[, variable_yaxis]
  dataset$variable_xaxis <- dataset[, variable_xaxis]
  
  # select age at time zero for which the plot should be created:
  # either median age (when age_plot not specified), or prespecified age
  if(!is.null(plot_age)){plot_age <- plot_age } else { plot_age = dataset$age_at_event %>% median()}
  
  # extract fixed effects part of the model formula
  fixed_effects <- paste0(formula(model))[3]
  model_formula_fixed <- paste0("~ ", fixed_effects)
  
  # define x-axis range for plotting (excluding extreme outliers)
  xaxis_lower <- quantile(dataset[!is.na(dataset$variable_yaxis), "variable_xaxis"], 0.005, na.rm = T) %>% as.numeric() %>% round(digits = 0)
  xaxis_upper <- 0 %>% as.numeric() %>% round(digits = 0)
  
  # create df_new  ----------------------------------------------------------------------------------------------------------
  df_new <- with(dataset, expand.grid(variable_xaxis = seq(xaxis_lower, xaxis_upper, length.out = 250),
                                      sex = levels(sex),
                                      event_status = levels(event_status),
                                      med_antihypertensive = levels(med_antihypertensive),
                                      med_lipidlowering = levels(med_lipidlowering),
                                      med_antidiabetic = levels(med_antidiabetic),
                                      age_at_event = plot_age,
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
    df_new <- with(dataset, expand.grid(variable_xaxis = seq(xaxis_lower, xaxis_upper, length.out = 250),
                                        sex = levels(sex),
                                        event_status = levels(event_status),
                                        med_antihypertensive = levels(med_antihypertensive),
                                        med_lipidlowering = levels(med_lipidlowering),
                                        med_antidiabetic = levels(med_antidiabetic),
                                        age_at_event = plot_age,
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
    mcoefs <- marginal_coefs(model, std_errors = TRUE, seed = 1)
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

trajectory_plot_tte <- function(df_new,
                                plot_antihypertensive = "No", plot_lipidlowering = "No",  plot_antidiabetic = "No",
                                plot_sex = "Female",
                                colors = c("#DD3C51", "#1F6683", "#313657"),
                                limits_xaxis = c(-35, 0), limits_yaxis = NULL,
                                plot_title = "", label_xaxis = "Time to event (years)", label_yaxis = "",
                                legend_labels = c("Individuals with primary outcome", "Survivors without primary outcome", "Individuals with non-cardiovascular death"),
                                y.ticks_multiplier = 2.5, y.at_multiplier = 5, legend_positions = c(1.01, 0.99, 0.97),
                                linetypes = c("solid", "solid", "solid"),
                                show_yaxis = TRUE, show_title = TRUE){
  
  # determine limits of xaxis and yaxis
  if(is.null(limits_xaxis)){
    xlim_lower <- round(min(df_new$time_to_event, na.rm = T),1)
    xlim_upper <- round(max(df_new$time_to_event, na.rm = T),1)
  } else {
    xlim_lower <- limits_xaxis[1]
    xlim_upper <- limits_xaxis[2]
  }
  
  # select 'subjects' that should be plotted
  df_plot <- subset(df_new, sex == plot_sex &
                      med_antihypertensive == plot_antihypertensive & med_lipidlowering == plot_lipidlowering & med_antidiabetic == plot_antidiabetic &
                      time_to_event >= xlim_lower & time_to_event <= xlim_upper)

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
  
  # split dataset for outcome
  df_survivors <- subset(df_plot, event_status == "Survivors without primary outcome")
  df_ncvd_death <- subset(df_plot, event_status == "Individuals with non-cardiovascular death")
  df_composite <- subset(df_plot, event_status == "Individuals with primary outcome")
  
  # plot  ----------------------------------------------------------------
  
  # plot (ggplot)
  plot <- ggplot() +
    # event_1
    geom_ribbon(data = df_ncvd_death, aes(x=time_to_event, ymin=lowerlimit, ymax=upperlimit), alpha=0.1, fill=colors[3]) +
    geom_line(data=df_ncvd_death, aes(x=time_to_event, y=pred), color=colors[3], linewidth = 1.25, linetype = linetypes[3]) +
    annotate("text", 
             label=legend_labels[3], 
             x=xlim_upper*0.45, 
             y=ylim_upper*legend_positions[3],
             col=colors[3], 
             hjust="right", 
             size = 5) +
    # event_2
    geom_ribbon(data = df_survivors, aes(x=time_to_event, ymin=lowerlimit, ymax=upperlimit), alpha=0.1, fill=colors[2]) +
    geom_line(data=df_survivors, aes(x=time_to_event, y=pred), color=colors[2], linewidth = 1.25, linetype = linetypes[2]) +
    annotate("text", 
             label=legend_labels[2], 
             x=xlim_upper*0.45, 
             y=ylim_upper*legend_positions[2],
             col=colors[2], 
             hjust="right", 
             size = 5) +
    # event_3
    geom_ribbon(data = df_composite, aes(x=time_to_event, ymin=lowerlimit, ymax=upperlimit), alpha=0.1, fill=colors[1]) +
    geom_line(data=df_composite, aes(x=time_to_event, y=pred), color=colors[1], linewidth = 1.25, linetype = linetypes[1]) +
    annotate("text", 
             label=legend_labels[1], 
             x=xlim_upper*0.45, 
             y=ylim_upper*legend_positions[1],
             col=colors[1], 
             hjust="right", 
             size = 5) +
    # titles and grid
    labs(x=label_xaxis, y=make_label(label_yaxis)) +
    ggtitle(plot_title) +
    theme(panel.background = element_rect(fill = "white", linewidth =  2, linetype = "solid"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey"), 
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size=16),
          axis.title.y = if(show_yaxis) element_text(size = 16) else element_text(color = NA),
          axis.text.y = if(show_yaxis) element_text(size = 16) else element_text(color = NA),
          axis.text.x = element_text(size=16),
          plot.title = if(show_title) element_text(size = 16, face = "bold", hjust = 0.5) else element_blank()) +
    scale_x_continuous(minor_breaks = x.ticks, breaks = c(x.at), expand=expansion(mult=c(0.05,0.05)))  +
    scale_y_continuous(minor_breaks = y.ticks, breaks = c(y.at), expand=expansion(mult=c(0.15,0.15))) +
    coord_cartesian(ylim=c(ylim_lower, ylim_upper), xlim=c(xlim_lower, xlim_upper))
  
  return(plot)
  
}

