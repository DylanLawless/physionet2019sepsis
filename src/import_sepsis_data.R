
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(forcats)
library(patchwork)  # For combining plots

# 	## Overview
# 	
# 	Data used in the competition is sourced from ICU patients in three separate hospital systems. Data from two hospital systems will be publicly available; however, one data set will be censored and used for scoring. The data for each patient will be contained within a single pipe-delimited text file. Each file will have the same header and each row will represent a single hour's worth of data. Available patient co-variates consist of Demographics, Vital Signs, and Laboratory values, which are defined in the tables below.
# 
# ## Time Points Definition
# 
# - **t suspicion:** Clinical suspicion of infection identified as the earlier timestamp of IV antibiotics and blood cultures within a specified duration. If antibiotics were given first, then the cultures must have been obtained within 24 hours. If cultures were obtained first, then antibiotic must have been subsequently ordered within 72 hours. Antibiotics must have been administered for at least 72 consecutive hours to be considered.
# - **t SOFA:** The occurrence of end organ damage as identified by a two-point deterioration in SOFA score within a 24-hour period.
# - **t sepsis:** The onset time of sepsis is the earlier of t suspicion and t SOFA as long as t SOFA occurs no more than 24 hours before or 12 hours after t suspicion; otherwise, the patient is not marked as a sepsis patient. Specifically, if \( t_{\text{suspicion}} - 24 \leq t_{\text{SOFA}} \leq t_{\text{suspicion}} + 12 \), then \( t_{\text{sepsis}} = \min(t_{\text{suspicion}}, t_{\text{SOFA}}) \).

# Import ----
# Set the paths for the raw data and the training set
path_raw_data <- "../data/physionet.org/files/challenge-2019/1.0.0/training/"
path_training_set <- "training_setA"
full_path <- paste0(path_raw_data, path_training_set)

# Get the list of .psv files and select the first 20 files
file_list <- list.files(path = full_path, pattern = "*.psv", full.names = TRUE)
selected_files <- file_list[1:500]

# Read the selected files and combine into a single dataframe
df_list <- lapply(selected_files, function(file) {
	read_data <- read.table(file, sep = "|", header = TRUE, fill = TRUE, quote = "")
	read_data$subject.id <- basename(file)
	return(read_data)
})

# Combine all data frames into one 'long' format dataframe
df <- do.call(rbind, df_list)

# Replace all "NaN" values with NA in the dataframe
df[sapply(df, is.nan)] <- NA

rm(file_list, full_path, path_raw_data, path_training_set)

# meta data ----
column_details <- data.frame(
	Category = c(rep("Vital Signs", 8), rep("Laboratory Values", 26), rep("Demographics", 6), "Outcome"),
	Column = c("HR", "O2Sat", "Temp", "SBP", "MAP", "DBP", "Resp", "EtCO2",
						 "BaseExcess", "HCO3", "FiO2", "pH", "PaCO2", "SaO2", "AST", "BUN", 
						 "Alkalinephos", "Calcium", "Chloride", "Creatinine", "Bilirubin_direct", 
						 "Glucose", "Lactate", "Magnesium", "Phosphate", "Potassium", 
						 "Bilirubin_total", "TroponinI", "Hct", "Hgb", "PTT", "WBC", "Fibrinogen", "Platelets",
						 "Age", "Gender", "Unit1", "Unit2", "HospAdmTime", "ICULOS",
						 "SepsisLabel"),
	Description = c("Heart rate (beats per minute)", "Pulse oximetry (%)", "Temperature (Deg C)", 
									"Systolic BP (mm Hg)", "Mean arterial pressure (mm Hg)", "Diastolic BP (mm Hg)", 
									"Respiration rate (breaths per minute)", "End tidal carbon dioxide (mm Hg)",
									"Measure of excess bicarbonate (mmol/L)", "Bicarbonate (mmol/L)", 
									"Fraction of inspired oxygen (%)", "N/A", 
									"Partial pressure of carbon dioxide from arterial blood (mm Hg)", 
									"Oxygen saturation from arterial blood (%)", "Aspartate transaminase (IU/L)", 
									"Blood urea nitrogen (mg/dL)", "Alkaline phosphatase (IU/L)", "Calcium (mg/dL)", 
									"Chloride (mmol/L)", "Creatinine (mg/dL)", "Bilirubin direct (mg/dL)", 
									"Serum glucose (mg/dL)", "Lactic acid (mg/dL)", "Magnesium (mmol/dL)", 
									"Phosphate (mg/dL)", "Potassium (mmol/L)", "Total bilirubin (mg/dL)", 
									"Troponin I (ng/mL)", "Hematocrit (%)", "Hemoglobin (g/dL)", 
									"Partial thromboplastin time (seconds)", "Leukocyte count (count*10^3/µL)", 
									"Fibrinogen (mg/dL)", "Platelets (count*10^3/µL)",
									"Years (100 for patients 90 or above)", "Female (0) or Male (1)", 
									"Administrative identifier for ICU unit (MICU)", "Administrative identifier for ICU unit (SICU)", 
									"Hours between hospital admit and ICU admit", "ICU length-of-stay (hours since ICU admit)",
									"For sepsis patients, SepsisLabel is 1 if t ≥ t_sepsis - 6 and 0 if t < t_sepsis - 6. For non-sepsis patients, SepsisLabel is 0.")
)

# View the dataframe
print(column_details)

