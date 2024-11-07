# Load required libraries
library(ggplot2)
library(survival)
library(KMsurv)
library(survminer)
library(moments)
library(dplyr)
library(car) # For VIF calculation
library(broom) # For tidy residuals

# Load the data
load("melanoma30.RData")
df_melanoma30 <- melanoma30

df_melanoma30$dead <- (as.numeric(df_melanoma30$dead)) - 1
df_melanoma30$epicell <- (as.numeric(df_melanoma30$epicell)) - 1
df_melanoma30$ulceration <- (as.numeric(df_melanoma30$ulceration)) - 1
df_melanoma30$sex <- (as.numeric(df_melanoma30$sex)) -1
df_melanoma30$invas2 <- (as.numeric(df_melanoma30$invas2)) - 1

# View the first few rows and structure of the dataset
head(df_melanoma30)
str(df_melanoma30)
summary(df_melanoma30)

# Boxplots for continuous variables
continuous_vars <- c("time", "thickness", "age", "logthick")
for (var in continuous_vars) {
  pdf(paste0("Billeder_duration/Boxplot_of_",var,".pdf"))
  ggplot(df_melanoma30, aes_string(x = var)) +
    geom_boxplot() +
    ggtitle(paste("Boxplot of", var)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(var) +
    ylab("Values") -> plot
  print(plot)
  dev.off()
}

# Histogram and density plot for continuous variables
for (var in continuous_vars) {
  pdf(paste0("Billeder_duration/Histogram_and_Density_of_",var,".pdf"))
  ggplot(df_melanoma30, aes_string(x = var)) +
    #  geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(aes(y = ..density.. * max(..count..)), color = "red") +
    ggtitle(paste("Histogram and Density of", var)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(var) -> plot
  print(plot)
  dev.off()
}

# Calculate skewness for continuous variables
sapply(df_melanoma30[continuous_vars], skewness, na.rm = TRUE)

# Table of frequencies for categorical variables
categorical_vars <- c("status", "dead", "ici", "epicell", "ulceration", "sex", "invas2")
lapply(df_melanoma30[categorical_vars], table)

# Bar plots for categorical variables
for (var in categorical_vars) {
  pdf(paste0("Billeder_duration/Bar_plot_of_",var,".pdf"))
  ggplot(df_melanoma30, aes_string(x = var)) +
    geom_bar(fill = "lightblue", color = "black") +
    ggtitle(paste("Bar plot of", var)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(var) -> plot
  print(plot)
  dev.off()
}





# Opgave 3: Kaplan-Meier Survival Analysis

# Create thickness categories
df_melanoma30$thickness_cat <- as.factor(cut(df_melanoma30$thickness,
                                             breaks = quantile(df_melanoma30$thickness, probs = c(0, 1/5, 2/5, 3/5, 4/5, 1), na.rm = TRUE),
                                             labels = c("Cat 1", "Cat 2", "Cat 3", "Cat 4", "Cat 5")))

# Create a survival object
df_melanoma30$dead <- as.numeric(df_melanoma30$dead)
survival_object <- Surv(df_melanoma30$time, df_melanoma30$dead)

# Fit Kaplan-Meier curves for thickness categories
fit <- survfit(survival_object ~ thickness_cat, data = df_melanoma30)

# Plot Kaplan-Meier survival curves
pdf("Billeder_duration/Kaplan-Meier_Survival_Curves_by_Tumor_Thickness_Categories.pdf")
ggsurvplot(fit, data = df_melanoma30, 
           #pval = TRUE, # Adds p-value from log-rank test
           #conf.int = FALSE, # Adds confidence intervals
           #risk.table = TRUE, # Shows risk table
           ggtheme = theme_minimal(), 
           palette = "Dark2",
           title = "Kaplan-Meier Survival Curves by Tumor Thickness Categories")
dev.off()
# Perform the log-rank test
log_rank_test <- survdiff(survival_object ~ thickness_cat, data = df_melanoma30)
log_rank_test





# Opgave 4: Cox Proportional Hazards Model

# Fit the Cox proportional hazards model
cox_model <- coxph(Surv(time, dead) ~ age + thickness + epicell + ici + ulceration + sex + invas2, data = df_melanoma30)
summary(cox_model)

# Calculate and add martingale residuals to the dataset
martingale_residuals <- residuals(cox_model, type = "martingale")
df_melanoma30$martingale_residuals <- martingale_residuals

# Plot martingale residuals against 'age'
pdf("Billeder_duration/Martingale_Residuals_vs_Age.pdf")
ggplot(df_melanoma30, aes(x = age, y = martingale_residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "blue") +
  ggtitle("Martingale Residuals vs Age") +
  xlab("Age") +
  ylab("Martingale Residuals")
dev.off()
# Plot martingale residuals against 'thickness'
pdf("Billeder_duration/Martingale_Residuals_vs_Tumor_Thickness.pdf")
ggplot(df_melanoma30, aes(x = thickness, y = martingale_residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "blue") +
  ggtitle("Martingale Residuals vs Tumor Thickness") +
  xlab("Tumor Thickness") +
  ylab("Martingale Residuals")
dev.off()


# Opgave 5 ---

# Assume you have already fitted your Cox model
# Example: cox_model <- coxph(Surv(time, dead) ~ age + thickness + epicell + ici + ulceration + sex + invas2, data = df_melanoma30)

## Check Proportional Hazards Assumption
ph_test <- cox.zph(cox_model)
print(ph_test)
jpeg("Billeder_duration/Schoenfeld_Residuals.jpeg")
# Use ggcoxzph to create a combined plot
ggcoxzph(ph_test)
# Close the device
dev.off()



## Examine Residuals
# Martingale Residuals
martingale_residuals <- residuals(cox_model, type = "martingale")
df_melanoma30$martingale_residuals <- martingale_residuals

# Plot martingale residuals against fitted values
pdf("Billeder_duration/martingale_residuals_against_fitted_values.pdf",width = 25,height = 12)
ggplot(df_melanoma30, aes(x = fitted(cox_model), y = martingale_residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "blue") +
  ggtitle("Martingale Residuals vs Fitted Values") +
  xlab("Fitted Values") +
  ylab("Martingale Residuals") +
  theme_minimal()
dev.off()
## Deviance Residuals
deviance_residuals <- residuals(cox_model, type = "deviance")
df_melanoma30$deviance_residuals <- deviance_residuals

# Plot deviance residuals
pdf("Billeder_duration/deviance_residuals.pdf")
ggplot(df_melanoma30, aes(x = fitted(cox_model), y = deviance_residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Deviance Residuals vs Fitted Values") +
  xlab("Fitted Values") +
  ylab("Deviance Residuals") +
  theme_minimal()
dev.off()

# Add Deviance Residuals to the dataframe
df_melanoma30$deviance_residuals <- deviance_residuals

# Plot Histogram of Deviance Residuals
devres=residuals(fitall,type="deviance")
plot(df_melanoma30$age,devres)
pdf("Billeder_duration/deviance_residuals_Histogram.pdf")
hist(devres)
dev.off()
boxplot(devres~df_melanoma30$dead)#many censored observations screw up deviance residuals.
mean(df_melanoma30$dead)
## Check for Influential Observations using DFBetas
# Calculate DFBetas for each observation in the Cox model
dfbetas_values <- residuals(cox_model, type = "dfbeta")

# Print the first few rows of the DFBetas to examine
head(dfbetas_values)

# Convert DFBetas matrix to a data frame for easier plotting
dfbetas_df <- as.data.frame(dfbetas_values)
names(dfbetas_df) <- names(coef(cox_model))  # Label columns with predictor names
dfbetas_df$Observation <- 1:nrow(dfbetas_df)  # Add observation index

# Plot DFBetas for each predictor
for (predictor in names(coef(cox_model))) { 
  pdf(paste0("Billeder_duration/DFBetas_for_",predictor,".pdf"))
  ggplot(dfbetas_df, aes_string(x = "Observation", y = predictor)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = c(-2/(sqrt(nrow(df_melanoma30))), 2/(sqrt(nrow(df_melanoma30)))), color = "red", linetype = "dashed") +  # Threshold lines
    ggtitle(paste("DFBetas for Predictor:", predictor)) +
    xlab("Observation") +
    ylab(paste("DFBetas for", predictor)) +
    theme_minimal() -> plot
  print(plot)
  dev.off()
}


## Cox snell
melanomasurv=Surv(df_melanoma30$time,df_melanoma30$dead)
#model with all covariates. Cox-Snell residuals
fitall <- coxph(melanomasurv ~ thickness + epicell 
                + ici + ulceration + sex 
                + invas2 + age, data=df_melanoma30)
summary(fitall)
resall <- residuals(fitall, type = "martingale")
#Cox-Snell residualer beregnes lettest udfra martingal-residualerne.
mart <- resall
coxsnell <- df_melanoma30$dead - mart
coxsnellfit <- survfit(Surv(coxsnell,df_melanoma30$dead)~1)



pdf("Billeder_duration/log_Cox_Snell_residualer.pdf")
plot(log(coxsnellfit$time),log(-log(coxsnellfit$surv)))
abline(c(0,1))
dev.off()
pdf("Billeder_duration/Cox_Snell_residualer.pdf")
plot(coxsnellfit$time,-log(coxsnellfit$surv))
abline(c(0,1))
dev.off()
#looks pretty good.


# Opgave 6 ----
cox_model_log <- coxph(Surv(time, dead) ~ log(age) + logthick + epicell + ici + ulceration + sex + invas2, data = df_melanoma30)
summary(cox_model_log)
#the p value of log thickness is not under so we cant reject the null hypothesis
#This suggests that logthick (tumor thickness) may not have a statistically
# significant effect on survival in this model.


# Opgave 7 ----
# Define the individual's covariate values
new_patient <- data.frame(
  logthick = log(750),
  epicell = 1,         # Assuming 1 represents presence of epitheloid cells
  ici = 2,
  ulceration = 1,      # Assuming 1 represents presence of ulceration
  sex = 1,             # Assuming 1 represents male (check your dataset coding)
  invas2 = 0,    #  I-III
  age = log(57)
)

# Calculate the baseline survival function
baseline_surv <- survfit(cox_model_log)

# Predict the survival curve for the new individual over time
surv_prob <- summary(survfit(cox_model_log, newdata = new_patient))

# Extract the survival probability at 8 years
time_8_years <- 8 * 365.25  # 8 years in days (if your time scale is in days)
surv_at_8_years <- surv_prob$surv[surv_prob$time == time_8_years]

# Find the closest time to 8 years in surv_prob$time
closest_time_index <- which.min(abs(surv_prob$time - time_8_years))

# Get the survival probability at the closest time point
surv_at_8_years <- surv_prob$surv[closest_time_index]

# Display the estimated survival probability at the closest time to 8 years
cat("The man's survival probability is", round(surv_at_8_years * 100, 1), "%")

# Generate the survival curve for the individual
surv_curve <- survfit(cox_model_log, newdata = new_patient)

# Plot the survidata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAkCAYAAAD7PHgWAAABBklEQVR4Xu2XMQrCQBBFBQvR6wgJHsEDpHVjBDvvoBhbI3bWCkZbFUyhFrYiEat0WgmC6AVkdQqbIVmWZAOi82C64b+/bDWZDEEQP4phTLMaa9d003bTGMgu1psF7JVGNzuWPdzs18GDz443rgrIcndXbvW8g1axGfZKo7P2eBXc+WB74a3FGXtiA1kwzfnpqTF7hL3SwDfAaz+BqvjkwYADe6WhglQwJlQwKVQwKakVTGOoYNL5z4JxwBlUMEwqAu9SwTCpCLxLBcOkIvCusoKT9/WFQ6OkIvCukoJwt5rO0sehUVIReBem6ng+OLBXmnKjn4PbGM5PeKnqgXIlo5vHXoL4Nl4ZYqbbEGA7+wAAAABJRU5ErkJggg==val curve
pdf("Billeder_duration/survival_curve_man_57.pdf")
plot(surv_curve, xlab = "Time (days)", ylab = "Survival Probability", main = "Survival Probability over Time",
     conf.int = TRUE, col = "blue", lwd = 2)
dev.off()



#Opgave 8 ----

# Estimate the cumulative baseline hazard from the Cox model
base_haz <- basehaz(cox_model_log, centered = FALSE)

# Plot the cumulative baseline hazard over time
pdf("Billeder_duration/cumulative_baseline_hazard_over_time.pdf")
plot(base_haz$time, base_haz$hazard, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Cumulative Baseline Hazard",
     main = "Cumulative Baseline Hazard from Cox Model")
dev.off()
# Plot log-log of cumulative hazard for assessing Weibull fit
plot(log(base_haz$time), log(base_haz$hazard), type = "l", col = "blue", lwd = 2,
     xlab = "log(Time)", ylab = "log(Cumulative Baseline Hazard)",
     main = "Log-Log Plot of Cumulative Baseline Hazard")

# Fit an exponential model
exp_model <- survreg(Surv(time, dead) ~ log(age) + logthick + epicell + ici + ulceration + sex + invas2, 
                     data = df_melanoma30, dist = "exponential")
summary(exp_model)
# Fit a Weibull model
weibull_model <- survreg(Surv(time, dead) ~ log(age) + logthick + epicell + ici + ulceration + sex + invas2, 
                         data = df_melanoma30, dist = "weibull")
summary(weibull_model)


