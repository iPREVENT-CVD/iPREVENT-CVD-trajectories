
#################################################################################################
# General information 
#################################################################################################

# Author: Marie de Bakker
# Created date: 11/09/2024
# Edited  16/09/2026 David Yeung

# Project: iPREVENT-CVD

# Information: Mixed effect models for cardiovascular risk factors

#################################################################################################
# Packages
#################################################################################################

library(dplyr)
library(nlme)
library(splines)
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



#################################################################################################
# Systolic blood pressure 
#################################################################################################


#  Set limits for boundary knots
bk <- c(quantile(dataset[!is.na(dataset$sbp) & !is.na(dataset$med_antihypertensive), "time_to_event"], 0.025, na.rm = T) %>% as.numeric(),
        quantile(dataset[!is.na(dataset$sbp) & !is.na(dataset$med_antihypertensive), "time_to_event"], 0.975, na.rm = T) %>% as.numeric())



fit_lme <- lme(fixed = as.formula(paste0("log2(sbp) ~ ns(time_to_event, df = 2, Boundary.knots = c(", bk[1], ",", bk[2], "))*event_status*sex*age_at_event*med_antihypertensive + cohort_name")),
               data = dataset,
               random = ~ 1 | study_subject_id,
               na.action = na.omit, 
               method = "REML",
               control = lmeControl(opt='optim'))


save(fit_lme, 
     file = "sbp_lme_event.Rdata")


#################################################################################################
# Diabetes (binomial model example)
#################################################################################################

bk <- c(quantile(dataset[!is.na(dataset$diabetes), "time_to_event"], 0.025, na.rm = T) %>% as.numeric(), 
        quantile(dataset[!is.na(dataset$diabetes), "time_to_event"], 0.975, na.rm = T) %>% as.numeric())


fit_lme <- mixed_model(fixed = as.formula(paste0("diabetes ~ ns(time_to_event, df = 3, Boundary.knots = c(", bk[1], ",", bk[2], "))*age_at_event*sex*event_status + cohort_name1 + cohort_name2 + cohort_name3")),
                       data = dataset,
                       random = ~ 1 | study_subject_id,
                       na.action = na.omit,
                       family = binomial(),
                       iter_EM = 0,
                       max_coef_value = 75)



# create and save predictions (predictions occur here for mixed model as boundary knots variable is required), age set as median ------
df_predictions <- predictions_for_plot(model = fit_lme,
                                       dataset = dataset,
                                       variable_yaxis = "diabetes",
                                       variable_xaxis = "time_to_event",
                                       plot_age = 69)


save(df_predictions, file = "diabetes_prediction_event.Rdata")

