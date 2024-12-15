# https://mc-stan.org/rstanarm/articles/jm.html

# Lets try
# Estimating Joint Models for Longitudinal and Time-to-Event Data with rstanarm

# The authors of the best result seem to predic sepsis based on the bias the sepsis cases have certain measurements present, rather than the observation values themselves. 
# >"One may be interested primarily in the evolution of the clinical biomarker but may wish to account for what is known as informative dropout. If the value of future (unobserved) biomarker measurements are related to the occurrence of the terminating event, then those unobserved biomarker measurements will be “missing not at random”"

library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default())

#  multiple-row per patient longitudinal biomarker information, as shown in
head(pbcLong)

# single-row per patient survival information, as shown in
head(pbcSurv)
d_pbcLong <- pbcLong
d_pbcSurv <- pbcSurv

# datasets
help("rstanarm-datasets", package = "rstanarm")

# Univariate joint model (current value association structure) ----

library(rstanarm)
mod1 <- stan_jm(formulaLong = logBili ~ sex + trt + year + (year | id),
								dataLong = pbcLong,
								formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt,
								dataEvent = pbcSurv,
								time_var = "year",
								chains = 1, refresh = 200, seed = 12345)

# We have a variety of methods and post-estimation functions available for this class, including: print, summary, plot, fixef, ranef, coef, VarCorr, posterior_interval, update, and more.

print(mod1)

# The output tells us that for each one unit increase in an individual’s underlying level of log serum bilirubin, their estimated log hazard of death increases by 35.1% (equivalent to a 3.86-fold increase in the hazard of death). The mean absolute deviation (MAD) is provided as a more robust estimate of the standard deviation of the posterior distribution. In this case the MAD_SD for the association parameter is 0.225, indicating there is quite large uncertainty around the estimated association between log serum bilirubin and risk of death (recall this is a small dataset containing only 40 patients!).

# If we wanted some slightly more detailed output for each of the model parameters, as well as further details regarding the model estimation (for example computation time, number of longitudinal observations, number of individuals, type of baseline hazard, etc) we can instead use the summary method:

summary(mod1, probs = c(.025,.975))


as.data.frame(VarCorr(mod1))


# Univariate joint model (current value and current slope association structure) ----
# In the previous example we were fitting a shared parameter joint model which assumed that the log hazard of the event (in this case the log hazard of death) at time t was linearly related to the subject-specific expected value of the longitudinal marker (in this case the expected value of log serum bilirubin) also at time t. This is the default association structure, although it could be explicitly specified by setting the assoc = "etavalue" argument.

# However, let’s suppose we believe that the log hazard of death is actually related to both the current value of log serum bilirubin and the current rate of change in log serum bilirubin. To estimate this joint model we need to indicate that we want to also include the subject-specific slope (at time t) from the longitudinal submodel as part of the association structure. We do this by setting the assoc argument equal to a character vector c("etavalue", "etaslope") which indicates our desired association structure:
	

mod2 <- stan_jm(formulaLong = logBili ~ sex + trt + year + (year | id),
								dataLong = pbcLong,
								formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt,
								dataEvent = pbcSurv,
								assoc = c("etavalue", "etaslope"),
								time_var = "year",
								chains = 1, refresh = 2000, seed = 12345)


# Multivariate joint model (current value association structures) ----

# Suppose instead that we were interested in two repeatedly measured clinical biomarkers, log serum bilirubin and serum albumin, and their association with the risk of death. We may wish to model these two biomarkers, allowing for the correlation between them, and estimating their respective associations with the log hazard of death. We will fit a linear mixed effects submodel (identity link, normal distribution) for each biomarker with a patient-specific intercept and linear slope but no other covariates. In the event submodel we will include gender (sex) and treatment (trt) as baseline covariates. Each biomarker is assumed to be associated with the log hazard of death at time 𝑡  via it’s expected value at time 𝑡 (i.e. a current value association structure).


mod3 <- stan_jm(
	formulaLong = list(
		logBili ~ sex + trt + year + (year | id),
		albumin ~ sex + trt + year + (year | id)),
	formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt,
	dataLong = pbcLong, dataEvent = pbcSurv,
	time_var = "year",
	chains = 1, refresh = 500, seed = 12345)

print(mod3)

summary(mod3, pars = "assoc")

# Posterior predictions ---
# We can also access the range of post-estimation functions (described in the stan_jm and related help documentation; see for example help(posterior_traj) or help(posterior_survfit)).
# 
# Predicted individual-specific longitudinal trajectory for in-sample individuals
# Predicted individual-specific biomarker values can be obtained using either the posterior_traj or posterior_predict function. The posterior_traj is preferable, because it can be used to obtain the biomarker values at a series of evenly spaced time points between baseline and the individual’s event or censoring time by using the default interpolate = TRUE option. Whereas, the posterior_predict function only provides the predicted biomarker values at the observed time points, or the time points in the new data. Predicting the biomarker values at a series of evenly spaced time points can be convenient because they can be easily used for plotting the longitudinal trajectory. Moreover, by default the posterior_traj returns a data frame with variables corresponding to the individual ID, the time, the predicted mean biomarker value, the limits for the 95% credible interval (i.e. uncertainty interval for the predicted mean biomarker value), and limits for the 95% prediction interval (i.e. uncertainty interval for a predicted biomarker data point), where the level for the uncertainty intervals can be changed via the prob argument. Conversely, the posterior_predict function returns an 𝑆 by 𝑁 matrix of predictions where 𝑆 is the number of posterior draws and 𝑁 is the number of prediction time points (note that this return type can also be obtained for posterior_traj by specifying the argument return_matrix = TRUE).
# 
# As an example, let’s plot the predicted individual-specific longitudinal trajectories for each of the two biomarkers (log serum bilirubin and serum albumin) in the multivariate joint model estimated above. We will do this for three individuals (IDs 6, 7 and 8) who were included in the model estimation.
# 
# Here are the plots for log serum bilirubin:
	
p1 <- posterior_traj(mod3, m = 1, ids = 6:8)
pp1 <- plot(p1, plot_observed = TRUE)
pp1

# and here are the plots for serum albumin:
	
p2 <- posterior_traj(mod3, m = 2, ids = 6:8)
pp2 <- plot(p2, plot_observed = TRUE)
pp2

# The m argument specifies which biomarker we want to predict for (only relevant for a multivariate joint model). The ids argument is optional, and specifies a subset of individuals for whom we want to predict. In the plotting method, the plot_observed = TRUE specifies that we want to include the observed biomarker values in the plot of the longitudinal trajectory.
# 
# If we wanted to extrapolate the trajectory forward from the event or censoring time for each individual, then this can be easily achieved by specifying extrapolate = TRUE in the posterior_traj call. For example, here is the plot for log serum bilirubin with extrapolation:

p3 <- posterior_traj(mod3, m = 1, ids = 6:8, extrapolate = TRUE)
pp3 <- plot(p3, plot_observed = TRUE, vline = TRUE)
pp3

# and for serum albumin with extrapolation:

p4 <- posterior_traj(mod3, m = 2, ids = 6:8, extrapolate = TRUE)
pp4 <- plot(p4, plot_observed = TRUE, vline = TRUE)
pp4

# Here, we included the vline = TRUE which adds a vertical dashed line at the timing of the individual’s event or censoring time. The interpolation and extrapolation of the biomarker trajectory can be further controlled through the control argument to the posterior_traj function; for example, we could specify the number of time points at which to predict, the distance by which to extrapolate, and so on.


# Predicted individual-specific survival curves for in-sample individuals ----
# Predicted individual-specific survival probabilities and/or survival curves can be obtained using the posterior_survfit function. The function by default returns a data frame with the individual ID, the time, and the predicted survival probability (posterior mean and limits for the 95% credible interval). The uncertainty level for the credible interval can be changed via the prob argument. By default, individual-specific survival probabilities are calculated conditional on the individual’s last known survival time. When we are predicting survival probabilities for individuals that were used in the estimation of the model (i.e. in-sample individuals, where no new covariate data is provided), then the individual’s “last known survival time” will be their event or censoring time. (Note that if we wanted didn’t want to condition on the individual’s last known survival time, then we could specify condition = FALSE, but we probably wouldn’t want to do this unless we were calculating marginal or standardised survival probabilities, which are discussed later).
# 
# The default argument extrapolate = TRUE specifies that the individual-specific conditional survival probabilities will be calculated at evenly spaced time points between the individual’s last known survival time and the maximum follow up time that was observed in the estimation sample. The behaviour of the extrapolation can be further controlled via the control argument. If we were to specify extrapolate = FALSE then the survival probabilities would only be calculated at one time point, which could be specified in the times argument (or otherwise would default to the individual’s last known survival time).
# 
# As an example, let plot the predicted individual-specific conditional survival curve for the same three individual’s that were used in the previous example. The predicted survival curve will be obtained under the multivariate joint model estimated above.

p5 <- posterior_survfit(mod3, ids = 6:8)
pp5 <- plot(p5)
pp5

# 
# We could customize the plot further, for example, by using any of the ggplot2 functionality or using the additional arguments described in help(plot.survfit.stanjm).
# 
# 
# Combined plot of longitudinal trajectories and survival curves ----
# The package also provides a convenience plotting function, which combines plots of the individual-specific longitudinal trajectories, and the individual-specific survival function. We can demonstrate this by replotting the predictions for the three individuals in the previous example:

plot_stack_jm(yplot = list(pp3, pp4), survplot = pp5)

# 
# Here we can see the strong relationship between the underlying values of the biomarkers and mortality. Patient 8 who, relative to patients 6 and 7, has a higher underlying value for log serum bilirubin and a lower underlying value for serum albumin at the end of their follow up has a far worse predicted probability of survival.
# 

# Predicted individual-specific longitudinal trajectory and survival curve for out-of-sample individuals (i.e. dynamic predictions) ----
# Let us take an individual from our training data, in this case the individual with subject ID value 8. However, we will pretend this individual was not a member of our training data and rather that they are a new individual for whom we have obtained new biomarker measurements. Our goal is to obtain predictions for the longitudinal trajectory for this individual, and their conditional survival curve, given that we know they are conditional on their biomarker measurements we currently have available.

# First, let’s extract the data for subject 8 and then rename their subject ID value so that they appear to be an individual who was not included in our training dataset:

ndL <- pbcLong[pbcLong$id == 8, , drop = FALSE]
ndE <- pbcSurv[pbcSurv$id == 8, , drop = FALSE]
ndL$id <- paste0("new_patient")
ndE$id <- paste0("new_patient")

# Note that we have both the longitudinal data and event data for this new individual. We require data for both submodels because we are going to generate dynamic predictions that require drawing new individual-specific parameters (i.e. random effects) for this individual conditional on their observed data. That means we need to evaluate the likelihood for the full joint model and that requires both the longitudinal and event data (note however that the status indicator death will be ignored, since it is assumed that the individual we are predicting for is still alive at the time we wish to generate the predictions).

# Now we can pass this data to the posterior_traj function in the same way as for the in-sample individuals, except we will now specify the newdataLong and newdataEvent arguments. We will also specify the last_time argument so that the function knows which variable in the event data specifies the individual’s last known survival time (the default behaviour is to use the time of the last biomarker measurement).

# Our predictions for this new individual for the log serum bilirubin trajectory can be obtained using:

p6 <- posterior_traj(mod3, m = 1,
											 newdataLong = ndL,
											 newdataEvent = ndE,
											 last_time = "futimeYears")

pp6 <- plot(p6, plot_observed = TRUE, vline = TRUE)
pp6

# and for the serum albumin trajectory:

p7 <- posterior_traj(mod3, m = 2,
											 newdataLong = ndL,
											 newdataEvent = ndE,
											 last_time = "futimeYears")

pp7 <- plot(p7, plot_observed = TRUE, vline = TRUE)
pp7


# For the conditional survival probabilities we use similar information, provided to the posterior_survfit function:
p8 <- posterior_survfit(mod3,
													newdataLong = ndL,
													newdataEvent = ndE,
													last_time = "futimeYears")

pp8 <- plot(p8)
pp8

# We can then use the plot_stack_jm function, as we saw in a previous example, to stack the plots of the longitudinal trajectory and the conditional survival curve:

plot_stack_jm(yplot = list(pp6, pp7), survplot = pp8)


# Here we see that the predicted longitudinal trajectories and conditional survival curve for this individual, obtained using the dynamic predictions approach, are similar to the predictions we obtained when we used their individual-specific parameters from the original model estimation. This is because in both situations we are conditioning on the same outcome data.
# 
# Side note: We can even compare the estimated individual specific parameters obtained under the two approaches. For example, here is the posterior mean for the estimated individual-specific parameters for individual 8 from the fitted model:

c(ranef(mod3)[["Long1"]][["id"]][8,],
	ranef(mod3)[["Long2"]][["id"]][8,])

# and here is the mean of the draws for the individual-specific parameters for individual 8 under the dynamic predictions approach:
colMeans(attr(p6, "b_new"))




# sepsis data ----

#  multiple-row per patient longitudinal biomarker information, as shown in
head(pbcLong)
# id      age sex trt      year     logBili albumin platelet


# single-row per patient survival information, as shown in
head(pbcSurv)
# id      age sex trt futimeYears status death


head(pbcLong)

# balance ----

# Count the number of sepsis cases
num_sepsis_cases <- sum(df$SepsisLabel == 1)

# Sample an equal number of non-sepsis cases
non_sepsis_sample <- df %>%
	filter(SepsisLabel == 0) %>%
	sample_n(num_sepsis_cases)

# Combine with sepsis cases
balanced_df <- df %>%
	filter(SepsisLabel == 1) %>%
	bind_rows(non_sepsis_sample)

# Check the balance
table(balanced_df$SepsisLabel)

# covert data to match package ----
df <- balanced_df 

# fill forward NAs
df_filled <- df %>%
	arrange(subject.id, ICULOS) %>%  # Arrange by subject.id and ICULOS
	group_by(subject.id) %>%  # Group by subject.id
	tidyr::fill(everything(), .direction = "down") %>%  # Fill NA with the previous non-NA value
	ungroup()  # Ungroup to prevent grouping affecting further operations

# Using dplyr to calculate and display the percentage of Non-NA values for each column
perc_non_na <- df %>%
	summarise(across(everything(), ~ 100 * sum(!is.na(.)) / nrow(df), .names = "perc_non_na_{.col}")) %>% round()

# Round the results
perc_non_na <- round(perc_non_na)

# Convert the data frame to a long format for easy filtering
library(tidyr)

perc_non_na_long <- perc_non_na %>%
	pivot_longer(cols = everything(), names_to = "column", values_to = "percentage")

# Filter to keep only columns with >= 80% non-NA values
columns_to_keep <- perc_non_na_long %>%
	filter(percentage >= 80) %>%
	pull(column)

# Extract just the names without the prefix for further processing if necessary
clean_column_names <- sub("perc_non_na_", "", columns_to_keep)

rm(perc_non_na, perc_non_na_long, columns_to_keep, clean_column_names)

names(df)

ggplot(df, aes(x = MAP)) +  # Use aes_string for dynamic variable names
	geom_histogram(bins = 30, fill = "#0e3b5c", color = "black")

ggplot(df, aes(x = Resp)) +  # Use aes_string for dynamic variable names
	geom_histogram(bins = 30, fill = "#0e3b5c", color = "black")


df <- df %>% select(
subject.id,
Age,
Gender,
SepsisLabel,
ICULOS,
MAP, 
O2Sat,
Resp)

# Extract unique and sorted subject.ids
unique_ids <- df %>%
	select(subject.id) %>%
	distinct() %>%
	arrange(subject.id) %>%
	mutate(id = row_number())  # Adds a sequential numeric ID

# View the mapping
print(unique_ids)


# Join the numeric IDs back to the original dataframe
df <- df %>%
	left_join(unique_ids, by = "subject.id")

# Check the first few rows of the updated dataframe
head(df)

df$age <- df$Age
df$sex <- df$Gender
df$year <- df$ICULOS
df$death <- df$SepsisLabel
df$logBili <- df$MAP
df$albumin <- df$Resp
df$platelet <- df$O2Sat
# df$trt <- 0

df <- df %>% select(
	-Age, 
	-subject.id,
	-Gender,
	-ICULOS,
	-SepsisLabel,
	# -Bilirubin_total,
	-Resp,
	# -Platelets
	)


names(df)

df <- df %>% na.omit()


# keep only data until a patient "death" (sepsis occurs)
df <- df %>%
	arrange(id, year) %>%  # Ensure data is ordered by id and year
	group_by(id) %>%  # Group by patient id
	mutate(cum_death = cumsum(death)) %>%  # Cumulative sum of death occurrences
	filter(cum_death < 2) %>%  # Keep only rows before and including the first death
	select(-cum_death)  # Remove the cumulative sum column



# create df_pbcSurv ----

# Sort by year.
# Keep the row with the first occurrence of death == 1.
# If there is no occurrence where death == 1, keep the row with the maximum year where death == 0.
# # Check the results
# head(df_pbcSurv)

df <- df %>%
	arrange(id, year) 

df_pbcSurv <- df %>%
	group_by(id) %>%
	mutate(last_year = max(year),  # Calculate the max year for each id
				 first_death_year = ifelse(any(death == 1), min(year[death == 1]), NA_real_)) %>%  # Calculate the first death year if exists
	filter((death == 1 & year == first_death_year) |  # Filter to keep the first death year
				 	(death == 0 & year == last_year & is.na(first_death_year))) %>%  # Or keep the last year if no death == 1
	select(-last_year, -first_death_year) %>%  # Clean up helper columns
	ungroup()  # Ensure no residual grouping affects further operations

head(df_pbcSurv)
nrow(df_pbcSurv)

df_pbcSurv$futimeYears <- df_pbcSurv$year
df_pbcSurv <- df_pbcSurv %>% select(-year)



# create df_pbcLong -----



#  multiple-row per patient longitudinal biomarker information, as shown in
head(pbcLong)
# id, age, sex, trt, year, logBili, albumin, platelet

df_pbcLong <- df %>% select(
	id, age, sex, year, logBili, albumin, platelet ) %>% distinct()
names(df )

# single-row per patient survival information, as shown in
head(pbcSurv)

# df_pbcSurv <- df %>% select(
# 	id, age, sex, trt, year, death ) %>% distinct()

# get the futimeYears where death occurred 
# # Assuming 'year' is your time variable and 'death' is the status indicator
# results <- df %>%
# 	group_by(id) %>%
# 	arrange(year) %>%  # Ensure data is ordered by year within each group
# 	mutate(change = death == 1 & lag(death, default = 0) == 0) %>%
# 	filter(change) %>%  # Filter to rows where a change from 0 to 1 occurs
# 	summarise(futimeYears = min(year))  # Get the minimum year of such change
# 
# library(dplyr)
# 
# # Calculate the maximum year where death is still 0 for each patient
# results <- df %>%
# 	group_by(id) %>%
# 	arrange(id, year) %>%  # Ensure data is ordered by year within each group
# 	filter(death == 0) %>%  # Keep only years where death is 0
# 	summarise(max_year_no_death = max(year, na.rm = TRUE))  # Get the maximum year where death is 0
# 
# # if death already occurred before ICU, then they would have NA, so lets find their min year
# min_death_years <- df %>%
# 	filter(death == 1) %>%
# 	group_by(id) %>%
# 	summarise(min_year_death = min(year, na.rm = TRUE))
# 
# # # Optionally join back with the original dataset for further analysis
# df_with_transition_year <- df %>%
# 	left_join(results, by = "id")
# 
# # Merge original df adding the year of death; when it occurs in ICU or min year if it occurred before ICU
# df_with_transition_year <- df_with_transition_year %>%
# 	left_join(min_death_years, by = "id") %>%
# 	mutate(futimeYears = ifelse(is.na(max_year_no_death), min_year_death, max_year_no_death)) %>%
# 	select(-min_year_death, -max_year_no_death)  # Optionally drop the temporary min_year_death column after updating
# 
# names(df_with_transition_year)
# 
# # # now filter down to the uniq set
# # df_pbcSurv <- df_with_transition_year %>% select(
# # 	id, age, sex, trt, futimeYears, death ) %>% 
# # 	# filter(death == 1) %>% 
# # 	distinct()
# 
# df_pbcSurv <- df_with_transition_year %>%
# 	select(id, age, sex, futimeYears, death) %>%
# 	arrange(id, death, futimeYears) %>%  # Ensure that death == 1 comes first and earlier years come first when death == 1
# 	group_by(id) %>%
# 	slice(1) %>%  # Take the first row per group after sorting
# 	ungroup()
# 
df$id %>% unique() %>% length()
df_pbcSurv$id %>% unique() %>% length()
df_pbcLong$id %>% unique() %>% length()
# # results$id  %>% unique() %>% length()

pbcLong <- df_pbcLong
pbcSurv <- df_pbcSurv 


## rerun analysis ----


# Multivariate joint model (current value association structures) ----

# Suppose instead that we were interested in two repeatedly measured clinical biomarkers, log serum bilirubin and serum albumin, and their association with the risk of death. We may wish to model these two biomarkers, allowing for the correlation between them, and estimating their respective associations with the log hazard of death. We will fit a linear mixed effects submodel (identity link, normal distribution) for each biomarker with a patient-specific intercept and linear slope but no other covariates. In the event submodel we will include gender (sex) and treatment (trt) as baseline covariates. Each biomarker is assumed to be associated with the log hazard of death at time 𝑡  via it’s expected value at time 𝑡 (i.e. a current value association structure).

# Checking factor levels in the longitudinal data
sapply(pbcLong, function(x) if(is.factor(x)) levels(x))

# Checking factor levels in the survival data
sapply(pbcSurv, function(x) if(is.factor(x)) levels(x))

pbcLong %>% as_tibble() %>% head()

head(pbcLong)
head(pbcSurv)

# Filter pbcLong Based on pbcSurv
# You need to join the two datasets on id and then filter out rows in pbcLong where year exceeds futimeYears for any given id.


# Join survival data to longitudinal data to compare times
# pbcLong <- pbcLong %>%
# 	left_join(pbcSurv %>% select(id, futimeYears, death), by = "id") %>%
# 	filter(year <= futimeYears) %>%
# 	select(-futimeYears)  # Remove the futimeYears column if no longer needed after filtering

# Check the first few rows of pbcLong to ensure the filtering is correct
head(pbcLong)



# drop trt
mod3 <- stan_jm(
	formulaLong = list(
		logBili ~ sex + year + (year | id),
		albumin ~ sex + year + (year | id)),
	formulaEvent = survival::Surv(futimeYears, death) ~ sex,
	dataLong = pbcLong, dataEvent = pbcSurv,
	time_var = "year",
	chains = 1, refresh = 500, seed = 12345)


# "number of observations (=390) <= number of random effects (=546) for term (year | id); the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable"
# This error arises because there are more levels in your random effects grouping than there are observations available to estimate these parameters, making the model overparameterized and the estimates of the random effects unidentifiable.


print("THIS IS SLOW !! 230 seconds")

mod3 <- stan_jm(
	formulaLong = list(
		logBili ~ sex + year + (1 | id),  # Random intercepts for `id`
		albumin ~ sex + year + (1 | id)
	),
	formulaEvent = survival::Surv(futimeYears, death) ~ sex,
	dataLong = pbcLong, dataEvent = pbcSurv,
	time_var = "year",
	chains = 1, refresh = 500, seed = 12345
)


print(mod3)
class(df$id)

summary(mod3, pars = "assoc")

# Posterior predictions ---
# We can also access the range of post-estimation functions (described in the stan_jm and related help documentation; see for example help(posterior_traj) or help(posterior_survfit)).
# 
# Predicted individual-specific longitudinal trajectory for in-sample individuals
# Predicted individual-specific biomarker values can be obtained using either the posterior_traj or posterior_predict function. The posterior_traj is preferable, because it can be used to obtain the biomarker values at a series of evenly spaced time points between baseline and the individual’s event or censoring time by using the default interpolate = TRUE option. Whereas, the posterior_predict function only provides the predicted biomarker values at the observed time points, or the time points in the new data. Predicting the biomarker values at a series of evenly spaced time points can be convenient because they can be easily used for plotting the longitudinal trajectory. Moreover, by default the posterior_traj returns a data frame with variables corresponding to the individual ID, the time, the predicted mean biomarker value, the limits for the 95% credible interval (i.e. uncertainty interval for the predicted mean biomarker value), and limits for the 95% prediction interval (i.e. uncertainty interval for a predicted biomarker data point), where the level for the uncertainty intervals can be changed via the prob argument. Conversely, the posterior_predict function returns an 𝑆 by 𝑁 matrix of predictions where 𝑆 is the number of posterior draws and 𝑁 is the number of prediction time points (note that this return type can also be obtained for posterior_traj by specifying the argument return_matrix = TRUE).
# 
# As an example, let’s plot the predicted individual-specific longitudinal trajectories for each of the two biomarkers (log serum bilirubin and serum albumin) in the multivariate joint model estimated above. We will do this for three individuals (IDs 6, 7 and 8) who were included in the model estimation.
# 
# Here are the plots for log serum bilirubin:

p1 <- posterior_traj(mod3, m = 1, ids = 6:8)
pp1 <- plot(p1, plot_observed = TRUE)
pp1

# and here are the plots for serum albumin:

p2 <- posterior_traj(mod3, m = 2, ids = 6:8)
pp2 <- plot(p2, plot_observed = TRUE)
pp2

# The m argument specifies which biomarker we want to predict for (only relevant for a multivariate joint model). The ids argument is optional, and specifies a subset of individuals for whom we want to predict. In the plotting method, the plot_observed = TRUE specifies that we want to include the observed biomarker values in the plot of the longitudinal trajectory.
# 
# If we wanted to extrapolate the trajectory forward from the event or censoring time for each individual, then this can be easily achieved by specifying extrapolate = TRUE in the posterior_traj call. For example, here is the plot for log serum bilirubin with extrapolation:

p3 <- posterior_traj(mod3, m = 1, ids = 6:8, extrapolate = TRUE)
pp3 <- plot(p3, plot_observed = TRUE, vline = TRUE)
pp3

# and for serum albumin with extrapolation:

p4 <- posterior_traj(mod3, m = 2, ids = 6:8, extrapolate = TRUE)
pp4 <- plot(p4, plot_observed = TRUE, vline = TRUE)
pp4

# Here, we included the vline = TRUE which adds a vertical dashed line at the timing of the individual’s event or censoring time. The interpolation and extrapolation of the biomarker trajectory can be further controlled through the control argument to the posterior_traj function; for example, we could specify the number of time points at which to predict, the distance by which to extrapolate, and so on.


# Predicted individual-specific survival curves for in-sample individuals ----
# Predicted individual-specific survival probabilities and/or survival curves can be obtained using the posterior_survfit function. The function by default returns a data frame with the individual ID, the time, and the predicted survival probability (posterior mean and limits for the 95% credible interval). The uncertainty level for the credible interval can be changed via the prob argument. By default, individual-specific survival probabilities are calculated conditional on the individual’s last known survival time. When we are predicting survival probabilities for individuals that were used in the estimation of the model (i.e. in-sample individuals, where no new covariate data is provided), then the individual’s “last known survival time” will be their event or censoring time. (Note that if we wanted didn’t want to condition on the individual’s last known survival time, then we could specify condition = FALSE, but we probably wouldn’t want to do this unless we were calculating marginal or standardised survival probabilities, which are discussed later).
# 
# The default argument extrapolate = TRUE specifies that the individual-specific conditional survival probabilities will be calculated at evenly spaced time points between the individual’s last known survival time and the maximum follow up time that was observed in the estimation sample. The behaviour of the extrapolation can be further controlled via the control argument. If we were to specify extrapolate = FALSE then the survival probabilities would only be calculated at one time point, which could be specified in the times argument (or otherwise would default to the individual’s last known survival time).
# 
# As an example, let plot the predicted individual-specific conditional survival curve for the same three individual’s that were used in the previous example. The predicted survival curve will be obtained under the multivariate joint model estimated above.

p5 <- posterior_survfit(mod3, ids = 6:8)
pp5 <- plot(p5)
pp5


p5 <- posterior_survfit(mod3, ids = 1:8)
pp5 <- plot(p5)
pp5

c <- posterior_survfit(mod3)
ppc <- plot(c)

print("can we just print the most intersting results? ")

# 
# We could customize the plot further, for example, by using any of the ggplot2 functionality or using the additional arguments described in help(plot.survfit.stanjm).
# 
# 
# Combined plot of longitudinal trajectories and survival curves ----
# The package also provides a convenience plotting function, which combines plots of the individual-specific longitudinal trajectories, and the individual-specific survival function. We can demonstrate this by replotting the predictions for the three individuals in the previous example:

plot_stack_jm(yplot = list(pp3, pp4), survplot = pp5)

Vital Signs
MAP
Mean arterial pressure (mm Hg)


Vital Signs
Resp
Respiration rate (breaths per minute)


df$logBili <- df$MAP
df$albumin <- df$Resp
df$platelet <- df$O2Sat




# 
# Here we can see the strong relationship between the underlying values of the biomarkers and mortality. Patient 8 who, relative to patients 6 and 7, has a higher underlying value for log serum bilirubin and a lower underlying value for serum albumin at the end of their follow up has a far worse predicted probability of survival.
# 

# Predicted individual-specific longitudinal trajectory and survival curve for out-of-sample individuals (i.e. dynamic predictions) ----
# Let us take an individual from our training data, in this case the individual with subject ID value 8. However, we will pretend this individual was not a member of our training data and rather that they are a new individual for whom we have obtained new biomarker measurements. Our goal is to obtain predictions for the longitudinal trajectory for this individual, and their conditional survival curve, given that we know they are conditional on their biomarker measurements we currently have available.

# First, let’s extract the data for subject 8 and then rename their subject ID value so that they appear to be an individual who was not included in our training dataset:

ndL <- pbcLong[pbcLong$id == 8, , drop = FALSE]
ndE <- pbcSurv[pbcSurv$id == 8, , drop = FALSE]
ndL$id <- paste0("new_patient")
ndE$id <- paste0("new_patient")

# Note that we have both the longitudinal data and event data for this new individual. We require data for both submodels because we are going to generate dynamic predictions that require drawing new individual-specific parameters (i.e. random effects) for this individual conditional on their observed data. That means we need to evaluate the likelihood for the full joint model and that requires both the longitudinal and event data (note however that the status indicator death will be ignored, since it is assumed that the individual we are predicting for is still alive at the time we wish to generate the predictions).

# Now we can pass this data to the posterior_traj function in the same way as for the in-sample individuals, except we will now specify the newdataLong and newdataEvent arguments. We will also specify the last_time argument so that the function knows which variable in the event data specifies the individual’s last known survival time (the default behaviour is to use the time of the last biomarker measurement).

# Our predictions for this new individual for the log serum bilirubin trajectory can be obtained using:

p6 <- posterior_traj(mod3, m = 1,
										 newdataLong = ndL,
										 newdataEvent = ndE,
										 last_time = "futimeYears")

pp6 <- plot(p6, plot_observed = TRUE, vline = TRUE)
pp6

# and for the serum albumin trajectory:

p7 <- posterior_traj(mod3, m = 2,
										 newdataLong = ndL,
										 newdataEvent = ndE,
										 last_time = "futimeYears")

pp7 <- plot(p7, plot_observed = TRUE, vline = TRUE)
pp7


# For the conditional survival probabilities we use similar information, provided to the posterior_survfit function:
p8 <- posterior_survfit(mod3,
												newdataLong = ndL,
												newdataEvent = ndE,
												last_time = "futimeYears")

pp8 <- plot(p8)
pp8

# We can then use the plot_stack_jm function, as we saw in a previous example, to stack the plots of the longitudinal trajectory and the conditional survival curve:

plot_stack_jm(yplot = list(pp6, pp7), survplot = pp8)


# Here we see that the predicted longitudinal trajectories and conditional survival curve for this individual, obtained using the dynamic predictions approach, are similar to the predictions we obtained when we used their individual-specific parameters from the original model estimation. This is because in both situations we are conditioning on the same outcome data.
# 
# Side note: We can even compare the estimated individual specific parameters obtained under the two approaches. For example, here is the posterior mean for the estimated individual-specific parameters for individual 8 from the fitted model:

c(ranef(mod3)[["Long1"]][["id"]][8,],
	ranef(mod3)[["Long2"]][["id"]][8,])

# and here is the mean of the draws for the individual-specific parameters for individual 8 under the dynamic predictions approach:
colMeans(attr(p6, "b_new"))


