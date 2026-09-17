
#################################################################################################
# General information 
#################################################################################################

# Author: Marie de Bakker
# Created date: 11/09/2024
# Edited  16/09/2026 David Yeung

# Project: iPREVENT-CVD

# Information: Mixed effect models for cardiovascular risk factors over lifecourse

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
source("function_lifecourse.R")

#################################################################################################
# Prepare dataset
#################################################################################################


# Order dataset by subject id and age
df_long <- df_long[order(df_long$study_subject_id, df_long$age),]


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



fit_lme <- lme(fixed = as.formula(paste0("log2(sbp) ~ ns(age, df = 2, Boundary.knots = c(", bk[1], ",", bk[2], "))sex*med_antihypertensive + cohort_name")),
               data = dataset,
               random = ~ 1 | study_subject_id,
               na.action = na.omit, 
               method = "REML",
               control = lmeControl(opt='optim'))


save(fit_lme, 
     file = "sbp_lme_age.Rdata")


#################################################################################################
# Diabetes (binomial model example)
#################################################################################################

# center age to allow model convergence
dataset <- dataset %>% filter(!is.na(diabetes)) %>% as.data.frame()

dataset$age_centered <- scale(dataset$age, center = T, scale = F)
mean_age <- mean(dataset$age, na.rm = T)



# Fitting of model
fit_lme <- mixed_model(fixed = diabetes ~ age_centered*sex + cohort_name1 + cohort_name2 + cohort_name3,
                       data = dataset,
                       random = ~ 1 | study_subject_id,
                       na.action = na.omit,
                       family = binomial(),
                       iter_EM = 0,
                       max_coef_value = 30)


# create and save predictions ------
df_predictions <- predictions_for_plot(model = fit_lme,
                                       dataset = dataset,
                                       variable_yaxis = "diabetes",
                                       variable_xaxis = "age_centered")

save(df_predictions, file = "diabetes_prediction_age.Rdata")

