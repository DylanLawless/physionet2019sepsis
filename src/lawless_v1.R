# https://mc-stan.org/rstanarm/articles/jm.html

# Lets try
# Estimating Joint Models for Longitudinal and Time-to-Event Data with rstanarm

# The authors of the best result seem to predic sepsis based on the bias the sepsis cases have certain measurements present, rather than the observation values themselves. 
# >"One may be interested primarily in the evolution of the clinical biomarker but may wish to account for what is known as informative dropout. If the value of future (unobserved) biomarker measurements are related to the occurrence of the terminating event, then those unobserved biomarker measurements will be “missing not at random”"


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

## Summary statistics ----

# Assuming 'df' is your main dataframe and contains the variables
# Retrieve all column names except 'subject.id'
unique_names <- names(df)[!names(df) %in% "subject.id"]

# Function to create a histogram for each unique name
create_histogram <- function(name, data) {
	subset_data <- data[[name]]  # Access the column directly
	
	if (anyNA(subset_data)) {
		warning(paste("NA values in", name))
	}
	
	ggplot(data, aes_string(x = name)) +  # Use aes_string for dynamic variable names
		geom_histogram(bins = 30, fill = "#0e3b5c", color = "black") +
		labs(title = name, y = "Count", x = name) +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

# Split names into batches for manageable plotting
name_chunks <- split(unique_names, ceiling(seq_along(unique_names) / 12))

# Generate and print plots
# Generate plots but do not print immediately
plot_batches <- lapply(name_chunks, function(chunk) {
	plots <- lapply(chunk, create_histogram, data = df)
	combined_plot <- wrap_plots(plots, ncol = 3)
	return(combined_plot)  # Just return, do not print
})

# Function to print plots
print_plots <- function(plot_list) {
	lapply(plot_list, print)
}

# Call this whenever you want to print the plots
print_plots(plot_batches)

# plot_batches[1]
# plot_batches[2]
# plot_batches[3]
# plot_batches[4]

# explore ----
# Calculate the max ICULOS (length of stay) for each subject and the maximum SepsisLabel (assuming SepsisLabel is binary 0 or 1)
subject_order <- df %>%
	group_by(subject.id) %>%
	summarise(
		MaxICULOS = max(ICULOS),
		MaxSepsisLabel = max(SepsisLabel)
	) %>%
	arrange(desc(MaxSepsisLabel), desc(MaxICULOS)) %>%
	pull(subject.id)

# Use fct_reorder to adjust levels of subject.id based on the calculated order
df$subject.id <- factor(df$subject.id, levels = subject_order)

# Create the tile plot
df %>%
	ggplot(aes(x = ICULOS, y = subject.id, fill = factor(SepsisLabel))) + 
	geom_tile(color = "black") +  # Create tiles
	scale_fill_manual(values = c("0" = "steelblue", "1" = "red")) + 
	labs(x = "ICU Length of Stay (hours)", 
			 y = "Subject ID", 
			 fill = "Sepsis Label") 


# test ----
class(df$SepsisLabel)
df$SepsisLabel %>% unique()

# Logistic Regression Model
# model_logistic <- glm(SepsisLabel ~ ., data = df, family = binomial())
model_logistic <- glm(SepsisLabel ~ ICULOS, data = df, family = binomial())

# Summary of the model to check p-values and coefficients
summary(model_logistic)

# Predicting probabilities
df$predicted_probability <- predict(model_logistic, type = "response")

# Visualising the predicted probabilities
ggplot(df, aes(x = ICULOS, y = predicted_probability, color = factor(SepsisLabel))) +
	geom_line(stat = "identity") +
	labs(title = "Predicted Probability of Sepsis over Time", x = "ICU Length of Stay (hours)", y = "Predicted Probability")


# test ----
# install.packages("rstanarm")
library(rstanarm)

# Fit a Bayesian logistic regression model
model_bayes <- stan_glm(SepsisLabel ~ Age + ICULOS, family = binomial, data = df, 
												prior = normal(0, 2.5), prior_intercept = normal(0, 2.5), 
												chains = 4, iter = 1000, seed = 123)

# Summarize the model
print(summary(model_bayes))

# Let's explore each section of this summary to better understand what the model tells us about the relationship between `Age`, `ICULOS` (time in ICU), and `SepsisLabel` (the binary outcome indicating the presence or absence of sepsis).
# ### 1. **Model Info**:
# - **Function**: `stan_glm` indicates the function used to fit the generalized linear model using Stan, which is a platform for statistical modeling and high-performance statistical computation.
# - **Family**: `binomial [logit]` specifies the family and link function used in the model. Here, it’s a binomial family with a logit link, appropriate for binary outcome data.
# - **Formula**: `SepsisLabel ~ Age + ICULOS` shows the model formula where `SepsisLabel` is predicted by `Age` and `ICULOS`.
# - **Algorithm**: `sampling` points out that the model uses sampling methods (specifically MCMC) to estimate the posterior distributions of the parameters.
# - **Sample**: `2000` is the posterior sample size used for inference after convergence.
# - **Priors**: This tells you the priors set for the model parameters (not detailed here but can be viewed using `help('prior_summary')`).
# - **Observations**: `19424` indicates the number of observations used in the model.
# - **Predictors**: `3` including the intercept (constant term), Age, and ICULOS.
# 
# ### 2. **Estimates**:
# Shows the posterior mean, standard deviation (sd), and percentile values (10%, 50%, 90%) for the model coefficients.
# - **(Intercept)**, **Age**, **ICULOS**:
#   - **Mean**: The posterior mean of the coefficients.
#   - **SD**: Standard deviation of the posterior distributions.
#   - Percentiles (10%, 50%, 90%): These give a range of plausible values for the coefficients at different percentiles, with 50% being the median.
# 
# ### 3. **Fit Diagnostics**:
# - **mean_PPD**: Stands for mean of the posterior predictive distribution of the outcome variable, which helps in checking how well the model predictions match the data.
# 
# ### 4. **MCMC Diagnostics**:
# Provide insight into the quality of the MCMC estimation.
# - **mcse**: Monte Carlo standard error, which measures the accuracy of the Monte Carlo estimation (smaller is better).
# - **Rhat**: The potential scale reduction factor, which should be close to 1.0 for convergence (values close to 1 suggest that multiple chains are converging to a common distribution).
# - **n_eff**: Effective sample size, providing a rough measure of the number of independent samples in the posterior distribution.
# 
# ### Visualizing the Output
# To visualize this summary and provide clearer insights, you can plot the posterior distributions of the coefficients and the model fit diagnostics:

library(rstan)
library(bayesplot)
library(rstanarm)

stan_plot(model_bayes, pars = c("Age", "ICULOS"))

# Extract posterior samples correctly
posterior_samples <- as.data.frame(model_bayes$stanfit)

# Setting a color scheme for better visualization
color_scheme_set("brightblue")

# Trace plots for 'Age' and 'ICULOS'
mcmc_trace(posterior_samples, pars = c("Age", "ICULOS"))

# Additional plots can include:
# Density plots for visualizing the distribution of the samples
mcmc_dens(posterior_samples, pars = c("Age", "ICULOS"))


# Assuming 'posterior_samples' already contains only the necessary columns for Age and ICULOS.
# If not, adjust the data frame accordingly:
posterior_long <- reshape2::melt(posterior_samples[, c("Age", "ICULOS")], variable.name = "Coefficient", value.name = "Value")

# Plotting
ggplot(posterior_long, aes(x = Value, fill = Coefficient)) +
	geom_density(alpha = 0.5) +
	scale_fill_manual(values = c("Age" = "blue", "ICULOS" = "red")) +
	labs(title = "Posterior Densities of Age and ICULOS Coefficients",
			 x = "Coefficient Value",
			 y = "Density") 

# Understanding the Difference in Means:
# 	Mean of Posterior for Age: This represents the average effect of a one-unit increase in Age (e.g., one year) on the log odds of developing sepsis, accounting for the uncertainty in the data and prior beliefs (as encoded by the prior distribution).
# 
# Mean of Posterior for ICULOS: Similarly, this is the average effect of a one-unit increase in ICULOS (e.g., one hour spent in the ICU) on the log odds of developing sepsis.
# 
# Interpretation:
# 	Greater Mean Value: A predictor (e.g., ICULOS) with a higher mean in its posterior distribution typically indicates a stronger association with the outcome variable, compared to predictors with lower means. If ICULOS has a higher mean than Age, this could suggest that time spent in the ICU is more strongly associated with the likelihood of developing sepsis than the patient’s age.
# 
# Sign and Magnitude: The sign of the mean tells you the direction of the association (positive or negative), while the magnitude gives an idea of the strength of this association. A positive coefficient implies an increase in the log odds of the outcome with an increase in the predictor, while a negative coefficient implies a decrease.
# 
# Overlap and Separation in Distributions: The amount of overlap between the distributions of two parameters can indicate the relative certainty or variability in their effects. Less overlap generally means more distinct effects on the outcome.

# R-hat plot for checking convergence across multiple chains
mcmc_rhat(posterior_samples)
names(posterior_samples)
class(posterior_samples$`(Intercept)`)

# Autocorrelation plot to check the independence of samples
mcmc_acf(posterior_samples, pars = c("Age", "ICULOS"))

# Checking model summary for overall diagnostics
print(summary(model_bayes))

# Launching shinystan for an interactive model exploration
library(shinystan)
launch_shinystan(model_bayes)
