

#################################################################################################
# General information 
#################################################################################################

# Author: Marie de Bakker
# Date: 13/05/2025
# Project: iPREVENT-CVD

# Information: plot mixed effect models (aging)

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
# Load data into R
#################################################################################################

load("A:/Datasets/Merged datasets/2025.09.12 df_iprevent and df_long_trajectories project.RData")

# functions
source("I:/PROJECTS/iPREVENT-CVD/Scripts/002_trajectories project/Functions/trajectory_plot_age.R")

#################################################################################################
# Preparations
#################################################################################################


# Order data by subject id and age
dataset <- df_long[order(df_long$study_subject_id, df_long$age),]

# set contrasts
contrasts(df_long$cohort_name) <- contr.sum(levels(df_long$cohort_name))


# Dummy groups for mixed model contrast
dataset$cohort_name1 <- ifelse(dataset$cohort_abbreviation == "CHS", 1, 
                               ifelse(dataset$cohort_abbreviation == "Whitehall", -1, 0))

dataset$cohort_name2 <- ifelse(dataset$cohort_abbreviation == "HUNT", 1, 
                               ifelse(dataset$cohort_abbreviation == "Whitehall", -1, 0))

dataset$cohort_name3 <- ifelse(dataset$cohort_abbreviation == "MESA", 1, 
                               ifelse(dataset$cohort_abbreviation == "Whitehall", -1, 0))


#################################################################################################
# SBP
#################################################################################################

# load model
load("sbp_lme_age.Rdata")

# create plot
df_predictions <- predictions_for_plot(model = fit_lme,
                                       dataset = dataset, 
                                       variable_yaxis = "sbp",
                                       variable_xaxis = "age")

plot_sbp <- trajectory_plot_age(df_new = df_predictions,
                            colors = c("#DD3C51", "#313657"),
                            legend_labels = c("", ""), 
                            plot_title = "Systolic blood pressure", limits_yaxis = c(110, 150),
                            label_xaxis = "Age (years)", label_yaxis = "Systolic blood pressure (mmHg)",
                            y.at_multiplier = 10,
                            limits_xaxis = c(20, 90))



#################################################################################################
# Diabetes
#################################################################################################

# load model
load("diabetes_prediction_age.Rdata")

# create plot
plot_diabetes <- trajectory_plot_age(df_new = df_predictions,
                                    colors = c("#DD3C51", "#313657"),
                                    legend_labels = c("", ""), 
                                    plot_title = "Diabetes", limits_yaxis = c(0,50),
                                     label_xaxis = "Age (years)", label_yaxis = "Estimated probability ('%' yes)",
                                    y.at_multiplier = 10,
                                    limits_xaxis = c(20, 90))


print(plot_sbp)

print(plot_diabetes)



