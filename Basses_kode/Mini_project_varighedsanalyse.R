# Load required libraries
library(ggplot2)
library(survival)
library(KMsurv)
library(survminer)
library(moments)
#stor tissemand
# Load the data
load("melanoma30.RData")
df_melanoma30 <- melanoma30

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
           pval = TRUE, # Adds p-value from log-rank test
           conf.int = FALSE, # Adds confidence intervals
           risk.table = TRUE, # Shows risk table
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
pdf()
# Plot martingale residuals against 'thickness'
pdf("Billeder_duration/Martingale_Residuals_vs_Tumor_Thickness.pdf")
ggplot(df_melanoma30, aes(x = thickness, y = martingale_residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "blue") +
  ggtitle("Martingale Residuals vs Tumor Thickness") +
  xlab("Tumor Thickness") +
  ylab("Martingale Residuals")
dev.off()
# Fit models with log-transformed age and thickness
cox_model_log <- coxph(Surv(time, dead) ~ log(age) + logthick + epicell + ici + ulceration + sex + invas2, data = df_melanoma30)
summary(cox_model_log)

# Compare AIC of original model vs log-transformed model
AIC(cox_model, cox_model_log)

# Opgave 5: Additional Steps (To be defined as needed)
# You can continue with model assessment, diagnostics, or any additional tasks for opgave 5
  
#Cox-Snell residuals: overall check of fit

# Martingale residuals: assessment of functional form of covariate

# Deviance residuals: detection of outliers

# Score-process residual: check of proportional hazards for each covariate

# Detection of influential observations. DF betas










