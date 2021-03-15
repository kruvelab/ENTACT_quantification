#libraries--------------------------------
library(car)
library(caret)
library(caTools)
library(dplyr)
library(enviPat)
library(gbm)
library(gsubfn)
library(janitor)
library(Metrics)
library(OrgMassSpecR)
library(plotly)
library(randomForest)
library(rcdk)
library(readr)
library(rcdklibs)
library(rgl)
library(rJava)
library(rjson)
library(RRF)
library(tidyverse)
library(webchem)
require(mgcv)

source("compound_eluent.R")
source("R theme.R")

#----------------------------------------reading in data ---------------------------------------
#presaved relevant descriptors and regressor
descs <-  read_rds("ESIpos_model_descs_191116.rds")
regressor <- read_rds("ESIpos_model_191116.rds")

#----------------------------reading the data-------------------------------------

# read the data with concentrations
lcms_data <-  read_delim('inputfiles/PosMode_Modeling_Data_Final_Training.csv',
                         delim = ",",
                         col_names = TRUE,
                         trim_ws = TRUE)

data_concentrations <- lcms_data %>%
  mutate(SMILES = MS_READY_SMILES) %>%
  gather('Corrected_HighConc', 'Corrected_MidConc', 'Corrected_LowConc', key = 'Concentration level', value = 'Concentration')
#concentration is uM

data_concentrations <- data_concentrations %>%
  mutate(Level = case_when(
    `Concentration level` == 'Corrected_HighConc' ~ 'High',
    `Concentration level` == 'Corrected_MidConc' ~ 'Mid',
    `Concentration level` == 'Corrected_LowConc' ~ 'Low'
  ))

data_signals <- lcms_data %>%
  gather('mean_High_Intensity', 'mean_Mid_Intensity', 'mean_Low_Intensity', key = 'Signal level', value = 'Signal')

data_signals <- data_signals %>%
  mutate(Level = case_when(
    `Signal level` == 'mean_High_Intensity' ~ 'High',
    `Signal level` == 'mean_Mid_Intensity' ~ 'Mid',
    `Signal level` == 'mean_Low_Intensity' ~ 'Low'
  )) %>%
  select(-c("std_High_Intensity", "std_Mid_Intensity", "std_Low_Intensity"))

dataset <- data_concentrations %>%
  left_join(data_signals) %>%
  select(-c('Concentration level', 'Signal level', 'Corrected_HighConc', 'Corrected_MidConc', 'Corrected_LowConc', 'mean_High_Intensity', 'mean_Mid_Intensity', 'mean_Low_Intensity', "std_High_Intensity", "std_Mid_Intensity", "std_Low_Intensity", "Ionization_Mode", "DTXSID", "Preferred_Name", "MS_READY_DTXCID", "MS_READY_FORMULA", "MS_READY_MONOISOTOPIC_MASS", "MS_READY_SMILES", "Observed_Isomer", "Isomer_Number", "Isomer_Order"))

ggplot(data = dataset) +
  geom_line(mapping = aes (x = log10(Concentration), y = log10(Signal), color = SMILES)) +
  theme(legend.position = "none") +
  my_theme


#read the desriptors
Padel_data <-  read_delim('inputfiles/descs_pos.csv',
                          delim = ",",
                          col_names = TRUE,
                          trim_ws = TRUE)

Padel_data = Padel_data %>%
  group_by(SMILES) %>%
  mutate(MW = molecularmass(SMILES),
         IC = isotopedistribution(SMILES)) %>%
  ungroup()

#joining the LCMS data with Padel descriptors
dataset <- dataset %>%
  left_join(Padel_data)

#reading the gradient program for reversed phase LC
eluent_parameters <- read_delim('inputfiles/eluent.csv',
                                delim = ",",
                                col_names = TRUE)
#organic modifier muutujana javast
organic_modifier <- "MeOH"
#pH saame ühe muutujana javast
pH <- 7
#NH4 iooni olemasolu saame muutujana javast
NH4 <- 1

#Combine all data togetehr and chaneg the names to the algorithm ready form
dataset <- dataset %>%
  mutate(
    c_M = Concentration /10^6,
    RF = Signal/c_M,
    organic_modifier = organic_modifier,
    organic = organicpercentage(eluent_parameters,mean_RT),
    pH.aq. = pH,
    NH4 = NH4,
    viscosity =  viscosity(organic,organic_modifier),
    surface_tension = surfacetension(organic,organic_modifier),
    polarity_index = polarityindex(organic,organic_modifier)) %>%
  #though all columns are selected with "everything()" in the end, the order of the columns is changed
  dplyr::select(SMILES,RF,c_M,Signal,mean_RT,organic,viscosity,surface_tension,polarity_index,everything())


#predictions based on the model developed in august
prediction_set_model_pos <- dataset %>%
  na.omit() %>%
  mutate(logIE_pred = 0)
prediction <-  predict(regressor, newdata = prediction_set_model_pos, predict.all = TRUE)
prediction <- prediction$aggregate
prediction_set_model_pos <- prediction_set_model_pos %>%
  mutate(logIE_pred = prediction) %>%
  select(SMILES,logIE_pred, everything())


#------------------------------predicting concentrations -----------------------------------

#linear regression for the calibration compounds
lin_fit_logRF <- lm(log(RF, 10) ~ logIE_pred, data = prediction_set_model_pos)

reg_plot<- lm(logIE_pred ~ log(RF, 10), data = prediction_set_model_pos)


ggplot(data = prediction_set_model_pos) +
  geom_point(mapping = aes(x = RF, y = logIE_pred, color = mean_RT), size = 3, alpha = 1/5) +
  ylim(0, 5) +
  scale_x_log10(limits = c(1e11, 1e15)) +
  geom_abline(slope = reg_plot$coefficients[2], intercept = reg_plot$coefficients[1], color='steelblue', size = 1) +
  scale_shape(solid = FALSE) +
  labs(x = "Response (Signal/concentration)", y = "predicted logIE") +
  my_theme


ggsave('Responce_to_IE_pos.svg')


#predicting for compounds which have concentration data available
prediction_set_model_pos <- prediction_set_model_pos %>%
  mutate(logRF_pred = lin_fit_logRF$coefficients[2]*logIE_pred + lin_fit_logRF$coefficients[1]) %>%
  mutate(c_pred = Signal/(10^logRF_pred)) %>%
  mutate(c_pred_acc =   case_when (
    c_pred > c_M ~ c_pred/c_M,
    c_pred < c_M ~ c_M/c_pred) ) %>%
  select(SMILES, c_pred, c_M, c_pred_acc, logRF_pred, logIE_pred, RF, everything())

#mean prediction error
mean(prediction_set_model_pos$c_pred_acc)
#3.05

#plot the predicted concentration against actual concentration
ggplot(data = prediction_set_model_pos) +
  geom_point(mapping = aes(x =c_M, y = c_pred, color = mean_RT, size = 10), alpha = 1/5) +
  scale_x_log10(limits = c(1e-9, 1e-5)) +
  scale_y_log10(limits = c(1e-9, 1e-5)) +
  geom_abline(slope = 1, intercept = 0, color='steelblue', size = 1) +
  scale_shape(solid = FALSE) +
  labs(x = "Spiked concentration (M)", y = "Predicted concentration (M)") +
  my_theme

ggsave("c_pred_pos.svg")

ggplot(data = prediction_set_model_pos) +
  geom_point(mapping = aes(x =c_M, y = c_pred_acc, color = mean_RT, size = 10), alpha = 1/5) +
  scale_x_log10(limits = c((10^-8.8), (10^-5.5))) +
  scale_y_log10(limits = c(1, 1e2)) +
  geom_abline(slope = 0, intercept = 3, color='steelblue', size = 1) +
  scale_shape(solid = FALSE) +
  labs(x = "Spiked concentration (M)", y = "Error in predicted concentration (times)") +
  my_theme

ggsave("c_pred_error_pos.svg")

#all data
dataset_pred <- dataset %>%
  select(SMILES, descs, viscosity, surface_tension, polarity_index, pH.aq., NH4) %>%
  mutate(logIE_pred = 0) %>%
  na.omit()

prediction_all <-  predict(regressor, newdata = dataset_pred, predict.all = TRUE)
prediction_all <- prediction_all$aggregate

dataset_pred <- dataset_pred %>%
  mutate(logIE_pred = prediction_all) %>%
  group_by(SMILES) %>%
  summarise(logIE_pred = mean(logIE_pred))

dataset <- dataset %>%
  select(SMILES, c_M, Level, RF, Signal) %>%
  left_join(dataset_pred) %>%
  select(SMILES, logIE_pred, RF, c_M, Signal, Level) %>%
  mutate(logRF_pred = lin_fit_logRF$coefficients[2]*logIE_pred + lin_fit_logRF$coefficients[1]) %>%
  mutate(c_pred = Signal/(10^logRF_pred)) %>%
  mutate(c_pred_acc =   case_when (
    c_pred > c_M ~ c_pred/c_M,
    c_pred < c_M ~ c_M/c_pred) ) %>%
  select(SMILES, c_pred, c_M, c_pred_acc, logRF_pred, logIE_pred, RF, Signal, Level)

ggplot(data = dataset) +
  geom_point(mapping = aes(x = factor(Level, c("Low", "Mid", "High")), y = c_pred), color = "steelblue", size = 4, alpha = 1/20) +
  scale_y_log10() +
  my_theme

ggsave("c_pred_all_pos.svg")

write_delim(dataset, "ENTACT_pos_predictions.csv", delim = ",")

#------evaluation after revealing the real concntrations--------


data_with_concnentration = read_delim("inputfiles/All_PosData_for_Evaluation.csv",
                                      delim = ",",
                                      col_names = TRUE)
data_with_concnentration = data_with_concnentration %>%
  mutate(Spike_Conc = Spike_Conc/10^6) %>%
  select(SMILES, Level, Spike_Conc)


data_combined = dataset %>%
  left_join(data_with_concnentration)

ggplot(data = data_combined) +
  geom_point(mapping = aes(x = log10(c_pred),
                           y = log10(Spike_Conc))) +
  geom_abline(intercept = 0, slope = 1) +
  my_theme

#predicting for compounds which have concentration data available
data_combined <- data_combined %>%
  filter(is.na(c_M)) %>%
  mutate(c_pred_acc =   case_when (
    c_pred > Spike_Conc ~ c_pred/Spike_Conc,
    c_pred < Spike_Conc ~ Spike_Conc/c_pred) ) %>%
  select(SMILES, c_pred, c_pred_acc, logRF_pred, logIE_pred, RF, everything())

#mean prediction error
mean(data_combined$c_pred_acc %>% na.omit())
#3.27
