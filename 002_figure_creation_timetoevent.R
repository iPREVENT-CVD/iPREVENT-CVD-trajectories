
#################################################################################################
# General information 
#################################################################################################

# Author: Marie de Bakker
# Date: 13/05/2025
# Project: iPREVENT-CVD

# Information: Plotting effects from mixed effect models


#################################################################################################
# Packages
#################################################################################################

library(dplyr)
library(nlme)
library(splines)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(cowplot)

#################################################################################################
# Load dataset into R
#################################################################################################

# Load dataset
load("df_long.RData")

# Load functions
source("function_timetoevent.R")


#################################################################################################
# Prepare dataset
#################################################################################################


# Order dataset by subject id and time to event
df_long <- df_long[order(df_long$study_subject_id, df_long$time_to_event),]


# Remove rows with data after event
dataset_data <- df_long[df_long$time_to_event <= 0 & 
                          !is.na(df_long$event_status), ]

# set contrasts
contrasts(dataset$cohort_name) <- contr.sum(levels(dataset$cohort_name))


# Dummy groups for mixed model contrast
dataset$cohort_name1 <- ifelse(dataset$cohort_abbreviation == "CHS", 1, 
                               ifelse(dataset$cohort_abbreviation == "Whitehall", -1, 0))

dataset$cohort_name2 <- ifelse(dataset$cohort_abbreviation == "HUNT", 1, 
                               ifelse(dataset$cohort_abbreviation == "Whitehall", -1, 0))

dataset$cohort_name3 <- ifelse(dataset$cohort_abbreviation == "MESA", 1, 
                               ifelse(dataset$cohort_abbreviation == "Whitehall", -1, 0))



# same color palette in males and in females
female_colors <- c("#DD3C51", "#DD3C51", "#999999")
male_colors <- c("#313657", "#313657", "#999999") 
linetypes = c("solid", "longdash", "solid")

#################################################################################################
# Plot age, median age is used for predictions for plots
#################################################################################################

df_iprevent %>% summarize(median = median(age_at_event),
                          twentyfive = quantile(age_at_event, 0.25),
                          seventyfive = quantile(age_at_event, 0.75),
                          eighty = quantile(age_at_event, 0.80))


plot_age = # Cohort median age

#################################################################################################
# SBP
#################################################################################################

# load model
load("sbp_lme_event.Rdata")

# create plot
df_predictions <- predictions_for_plot(model = fit_lme,
                                       dataset = dataset,
                                       variable_yaxis = "sbp",
                                       variable_xaxis = "time_to_event",
                                       plot_age = plot_age)

plot_F_sbp <- trajectory_plot_tte(df_new = df_predictions,
                              plot_sex = "Female",
                              plot_title = "Systolic blood pressure", 
                              label_xaxis = "Time to event (years)", label_yaxis = "Systolic blood pressure (mmHg)",
                              limits_yaxis = c(115, 140), y.at_multiplier = 5, limits_xaxis = c(-30, 0), 
                              legend_positions = c(2, 2, 2),
                              colors = female_colors,
                              linetypes = linetypes,
                              show_yaxis = TRUE,
                              show_title = FALSE)

plot_M_sbp <- trajectory_plot_tte(df_new = df_predictions,
                              plot_sex = "Male",
                              plot_title = "Systolic blood pressure", 
                              label_xaxis = "Time to event (years)", label_yaxis = "Systolic blood pressure (mmHg)",
                              limits_yaxis = c(115, 140),y.at_multiplier = 5, limits_xaxis = c(-30, 0),
                              legend_positions =  c(2, 2, 2),
                              colors = male_colors,
                              linetypes = linetypes,
                              show_yaxis = FALSE,
                              show_title = FALSE)

#################################################################################################
# Diabetes
#################################################################################################

# load model
load("diabetes_prediction_event.Rdata")

plot_F_diabetes <- trajectory_plot_tte(df_new = df_predictions,
                                      plot_sex = "Female",
                                      plot_title = "Diabetes", 
                                      label_xaxis = "Time to event (years)", label_yaxis = "Diabetes ('%' yes)",
                                      limits_yaxis = c(0, 20), y.at_multiplier = 5, limits_xaxis = c(-30, 0), 
                                      legend_positions = c(2, 2, 2),
                                      colors = female_colors,
                                      linetypes = linetypes,
                                      show_yaxis = TRUE,
                                      show_title = FALSE)

plot_M_diabetes <- trajectory_plot_tte(df_new = df_predictions,
                                      plot_sex = "Male",
                                      plot_title = "Diabetes", 
                                      label_xaxis = "Time to event (years)", label_yaxis = "Diabetes ('%' yes)",
                                      limits_yaxis = c(0, 20), y.at_multiplier = 5, limits_xaxis = c(-30, 0),
                                      legend_positions =  c(2, 2, 2),
                                      colors = male_colors,
                                      linetypes = linetypes,
                                      show_yaxis = FALSE,
                                      show_title = FALSE)

#################################################################################################
# Combine plots for each risk factor and print
#################################################################################################

plot_sbp  <- ggdraw() + 
  draw_plot(plot_grid(plot_F_sbp, plot_M_sbp, ncol = 2, rel_widths = c(1, 1)), y = 0, height = 0.9) +
  draw_label("Systolic blood pressure", x = 0.5, y = 0.95, hjust = 0.5, fontface = "bold", size = 22)

plot_diabetes   <- ggdraw() + 
  draw_plot(plot_grid(plot_F_diabetes, plot_M_diabetes, ncol = 2, rel_widths = c(1, 1)), y = 0, height = 0.9) +
  draw_label("Diabetes", x = 0.5, y = 0.95, hjust = 0.5, fontface = "bold", size = 22)

print(plot_sbp)

print(plot_diabetes)

