load("C:/Users/jensp/Downloads/melanoma30.RData")


library(ggplot2)
library(survival)
library(KMsurv)
library(survminer)
library(moments)

df_melanoma30 <- melanoma30

# View the first few rows and structure
head(df_melanoma30)
str(df_melanoma30)
summary(df_melanoma30)


# Boxplots for continuous variables
continuous_vars <- c("time", "thickness", "age", "logthick")
for (var in continuous_vars) {
  ggplot(df_melanoma30, aes_string(x = var)) +
    geom_boxplot() +
    ggtitle(paste("Boxplot of", var)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(var) +
    ylab("Values") -> plot
  print(plot)
}

# Histogram and density plot for continuous variables
for (var in continuous_vars) {
  ggplot(df_melanoma30, aes_string(x = var)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(aes(y = ..density.. * max(..count..)), color = "red") +
    ggtitle(paste("Histogram and Density of", var)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(var) -> plot
  print(plot)
}


# Calculate skewness for continuous variables
sapply(df_melanoma30[continuous_vars], skewness, na.rm = TRUE)


# Table of frequencies for categorical variables
categorical_vars <- c("status", "dead", "ici", "epicell", "ulceration", "sex", "invas2")

# Frequency tables
lapply(df_melanoma30[categorical_vars], table)

# Bar plots for categorical variables
for (var in categorical_vars) {
  ggplot(df_melanoma30, aes_string(x = var)) +
    geom_bar(fill = "lightblue", color = "black") +
    ggtitle(paste("Bar plot of", var)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(var) -> plot
  print(plot)
}

#opgave 3

df_melanoma30$thickness_cat <-  as.factor(cut(df_melanoma30$thickness,
                                    breaks = quantile(df_melanoma30$thickness, probs = c(0, 1/5, 2/5, 3/5, 4/5, 1), na.rm = TRUE),
                                    labels = c("Cat 1", "Cat 2", "Cat 3", "Cat 4", "Cat 5")))

df_melanoma30$dead <- as.numeric(df_melanoma30$dead)
survival_object <- Surv(df_melanoma30$time, df_melanoma30$dead)

# Fit Kaplan-Meier survival curves for thickness categories
fit <- survfit(survival_object ~ thickness_cat, data = df_melanoma30)

# Plot the survival curves
ggsurvplot(fit, data = df_melanoma30, 
           pval = TRUE, # Adds p-value from log-rank test
           conf.int = TRUE, # Adds confidence intervals
           risk.table = TRUE, # Shows risk table
           ggtheme = theme_minimal(), 
           palette = "Dark2",
           title = "Kaplan-Meier Survival Curves by Tumor Thickness Categories")

#plot(fit, col = 1:3, main = "Kaplan-Meier Survival Curves by Tumor Thickness Categories")
#legend("topright", legend = levels(df_melanoma30$thickness_cat), col = 1:3, lty = 1)


# Perform the log-rank test independently
log_rank_test <- survdiff(survival_object ~ thickness_cat, data = df_melanoma30)

# Print the log-rank test result
log_rank_test

#opgave 4

# Fit the Cox proportional hazards model
cox_model <- coxph(Surv(time, dead) ~ age + thickness + epicell + ici 
                   + ulceration + sex + invas2, data = df_melanoma30)
summary(cox_model)


# Calculate martingale residuals
martingale_residuals <- residuals(cox_model, type = "martingale")

# Add martingale residuals to the dataframe
df_melanoma30$martingale_residuals <- martingale_residuals

# Plot martingale residuals for 'age'
ggplot(df_melanoma30, aes(x = age, y = martingale_residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "blue") +
  ggtitle("Martingale Residuals vs Age") +
  xlab("Age") +
  ylab("Martingale Residuals")

# Plot martingale residuals for 'thickness'
ggplot(df_melanoma30, aes(x = thickness, y = martingale_residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "blue") +
  ggtitle("Martingale Residuals vs Tumor Thickness") +
  xlab("Tumor Thickness") +
  ylab("Martingale Residuals")

# Fit a model with log-transformed age and thickness
cox_model_log <- coxph(Surv(time, dead) ~ log(age) + log(thickness) + epicell + ici 
                       + ulceration + sex + invas2, data = df_melanoma30)

cox_model_log <- coxph(Surv(time, dead) ~ log(age) + log(thickness), data = df_melanoma30)
cox_model <- coxph(Surv(time, dead) ~ age + thickness, data = df_melanoma30)
# Compare AIC of original model vs log-transformed model
AIC(cox_model, cox_model_log)



# Opgave 5















