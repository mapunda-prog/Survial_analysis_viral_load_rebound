#importing packages
library(readxl) # for reading excel files
library(tidyverse) # for data manipulation
library(survival) # for survival analysis
library(survminer) # for plotting survival curves
library(rio) # for importing data
library(janitor)  # for clean_names function
library(ggsurvfit) #creating survival plots
library(gtsummary)
library(glue)


# loading packages for data cleaning ----

pacman::p_load(
  rio,        # importing data  
  here,       # relative file pathways  
  janitor,    # data cleaning and tables
  lubridate,  # working with dates
  matchmaker, # dictionary-based cleaning
  epikit,     # age_categories() function
  tidyverse   # data management and visualization
)




#importing data set

Rebound_data <- import("data/combined_data2.csv") 

# explore dataset
names(Rebound_data)

# Select needed columns for analysis
Rebound_data_clean <- Rebound_data %>% 
  select("facility_council", "facility_name", 
                                              "patient_id", "sex", "age", "art_start_date",
                                              "marital_status", "residence", "residence_lvl", 
                                              "distance_HF", "stigma", "mental_illness", "alcohol_use",
                                              "cd4_test_date", "cd4_counts", "who_clinical_stage", 
                                              "art_adherence", "art_regimen", "arv_code", "arv_combination", 
                                              "hiv_care_appointment", "number_of_viral_load_test", "baseline_hvl_date", 
                                              "HVL_baseline", "baseline_hvl_copies_ml", "low_level_viremia_at_basline", 
                                              "duration_on_art_months", "art_dispensing_day", "hiv_tb_co_infection",
                                              "tb_preventive_therapy", "interruption_in_art_treatment", "client_category",
                                              "facility_ownership", "type_of_health_facility", "facility_patient_volume", 
                                              "facility_staffing", "art_refill_model", "Rebound_ovall", "DAYZ_REBOUND")
# create age categories ----
Rebound_data_clean <- Rebound_data_clean %>%
  mutate(age_cat = age_categories(age, 
                                  breakers = c(15, 25, 45,65, 85),
                                  ceiling = TRUE)) # 85 is ceiling, all above become NA

# get a summary plot for age categories
A <- plot(Rebound_data_clean$age_cat)


#inspect the data set
glimpse(Rebound_data_clean)

#changing the character variables into factor
Rebound_data_clean <- Rebound_data_clean %>% 
  mutate_if(is.character, as.factor)


# Creating the survival object

Surv(Rebound_data_clean$DAYZ_REBOUND, Rebound_data_clean$Rebound_ovall)[1:100]

km_Fisher <- survfit(Surv(DAYZ_REBOUND, Rebound_ovall) ~ 1, data = Rebound_data_clean)

# Create the survival curve plot
ggsurvplot(km_Fisher, data = Rebound_data_clean, risk.table = TRUE, 
           surv.median.line = "hv", censor = FALSE)

str(km_Fisher)

#ploting the survival function
km_Fisher2 <- ggsurvfit::survfit2(Surv(DAYZ_REBOUND, Rebound_ovall) ~ 1, data = Rebound_data_clean) %>%
  ggsurvfit() +
  labs(
    x = "Days to viral load rebound",
    y = "Overall rebound probability"
  ) + 
  add_confidence_interval() + add_risktable()

km_Fisher2

# Estimating x-year survival
# probability of surviving rebound beyond a certain number of years
summary(km_Fisher, times = 365.25) 



#Univariate Cox regression model

coxph(Surv(DAYZ_REBOUND, Rebound_ovall) ~ sex, data = Rebound_data_clean) %>%
  tbl_regression(exp = TRUE) 
coxph(Surv(DAYZ_REBOUND, Rebound_ovall) ~ age_cat, data = Rebound_data_clean) %>%
  tbl_regression(exp = TRUE) 
coxph(Surv(DAYZ_REBOUND, Rebound_ovall) ~ residence_lvl, data = Rebound_data_clean) %>%
  tbl_regression(exp = TRUE)

coxph(Surv(DAYZ_REBOUND, Rebound_ovall) ~ art_regimen, data = Rebound_data_clean) %>%
  tbl_regression(exp = TRUE)

coxph(Surv(DAYZ_REBOUND, Rebound_ovall) ~ hiv_care_appointment, data = Rebound_data_clean) %>%
  tbl_regression(exp = TRUE)

coxph(Surv(DAYZ_REBOUND, Rebound_ovall) ~ client_category, data = Rebound_data_clean) %>%
  tbl_regression(exp = TRUE)

coxph(Surv(DAYZ_REBOUND, Rebound_ovall) ~ art_refill_model, data = Rebound_data_clean) %>%
  tbl_regression(exp = TRUE)


# Export data for analysis in other tools
export(Rebound_data_clean, "Rebound.csv")

