# Thomas Leahy
# 23/01/2021
# Question 1


# Packages ----------------------------------------------------------------
if (!require("glmnet")) install.packages("glmnet")
if (!require("scales")) install.packages("scales")
if (!require("ggplot2")) install.packages("ggplot2")


# Data + format + transform -----------------------------------------------

setwd("C:/Users/thoma/OneDrive/Desktop/SH_assignment")
raw_data <- read.table("processed.cleveland.data", sep = ",", stringsAsFactors = F)
colnames(raw_data) <-c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",  "thalach",
                   "exang", "oldpeak", "slope", "ca","thal", "num")

#format missing data as NA
raw_data <- replace(raw_data, raw_data == "?", NA)

# regrouping of variables
raw_data$restecg[raw_data$restecg == 2] <- 1
raw_data$ca[raw_data$ca == "3.0"] <- "2.0"

# set as factor variables
cols_factor <- c("sex", "cp", "fbs", "restecg", "exang","slope", "thal", "ca")
raw_data[cols_factor] <- lapply(raw_data[cols_factor], factor)
raw_data$sex <- factor(raw_data$sex, levels = c(0,1), labels = c("Female", "Male"))

# remove missing data (<2%)
data_compcase <- na.omit(raw_data)

# convert outcome varible
data_compcase$outcome <- ifelse(data_compcase$num==0,0,1)
data_compcase <- subset.data.frame(data_compcase, select = -c(num))


# Data summary ------------------------------------------------------------

table(data_compcase$outcome, data_compcase$sex)
table(data_compcase$outcome, data_compcase$cp)
table(data_compcase$outcome, data_compcase$fbs)
table(data_compcase$outcome, data_compcase$restecg)
table(data_compcase$outcome, data_compcase$exang)
table(data_compcase$outcome, data_compcase$slope)
table(data_compcase$outcome, data_compcase$thal)
table(data_compcase$outcome, data_compcase$ca)
aggregate(data_compcase[,c("age")], list(data_compcase$outcome), mean)
aggregate(data_compcase[,c("age")], list(data_compcase$outcome), sd)
aggregate(data_compcase[,c("chol")], list(data_compcase$outcome), mean)
aggregate(data_compcase[,c("chol")], list(data_compcase$outcome), sd)
aggregate(data_compcase[,c("trestbps")], list(data_compcase$outcome), mean)
aggregate(data_compcase[,c("trestbps")], list(data_compcase$outcome), sd)
aggregate(data_compcase[,c("thalach")], list(data_compcase$outcome), mean)
aggregate(data_compcase[,c("thalach")], list(data_compcase$outcome), sd)
aggregate(data_compcase[,c("oldpeak")], list(data_compcase$outcome), mean)
aggregate(data_compcase[,c("oldpeak")], list(data_compcase$outcome), sd)


# Model -------------------------------------------------------------------

# Univariate analysis#

# fit logitic models
uni <- lapply(colnames(data_compcase)[1:13],
              
              function(var) {
                
                formula    <- as.formula(paste("outcome ~", var))
                logistic.fit <- glm(formula, data = data_compcase, family = binomial)
                
              })
# summarise ORs for fitted models
univar_res <- data.frame()
for (u in 1:13){
  for (est in 2:length(coef(uni[[u]]))){
    row <- exp(cbind(coef(uni[[u]]), confint(uni[[u]])))[est,]  
    name <- names(coef(uni[[u]]))[est]
    pval <- summary(uni[[u]])$coefficients[est,4]  
    add_row <- c(name, row,pval)
    univar_res <- rbind.data.frame(univar_res, add_row, stringsAsFactors = F)
  }
  
}
colnames(univar_res) <- c("Cov", "OR", "low95", "upp95", "pval")
cols_numeric <- c( "OR", "low95", "upp95")
univar_res[cols_numeric] <- lapply(univar_res[cols_numeric], as.numeric)

# Multicollinearity  #
cor_data <- data.matrix(subset.data.frame(data_compcase, select = -c(outcome)))

cors <- matrix(nrow = 13, ncol = 13)
for (i in 1:13){
  for (j in 1:13){
    test_cor <- cor.test(cor_data[,i], cor_data[,j])
    cors[i,j] <- as.numeric(test_cor$estimate)
  }
}
colnames(cors) <- colnames(cor_data)
rownames(cors) <- colnames(cor_data)

# Multivariate #
full_data <- subset.data.frame(data_compcase, select = -c(oldpeak))
full_model <- glm(outcome ~ ., data = full_data, family = binomial)
summary(full_model)
full_model_res <- as.data.frame(exp(cbind(coef(full_model), confint(full_model))))
full_model_res <- tail(full_model_res, -1)
full_model_res <- cbind.data.frame(rownames(full_model_res), full_model_res)
colnames(full_model_res) <- c("var", "or", "cilo","ciup" )

# LASSO #

x_vars <- model.matrix(outcome ~. , full_data)[,-1]
y_var <- full_data$outcome

# cv to selct lambda
cv_output <- cv.glmnet(x_vars, y_var,
                       alpha = 1, nfolds = 5)
# identifying best lamda
lambda_1se <- cv_output$lambda.1se
lambda_1se

# final lasso
lasso_best <- glmnet(x_vars, y_var, alpha = 1, lambda = lambda_1se)
coef(lasso_best)


# Plot --------------------------------------------------------------------

# Forest plot univariate
fp_univariate <- ggplot(data=univar_res, aes(x=Cov, y=OR, ymin=low95, ymax=upp95)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariate") + ylab("Odd ratio (95% CI)") +
  theme_bw() # use a white background 
fp_univariate  <- fp_univariate + scale_y_continuous(trans = log_trans(), breaks = c(0.2,0.5,1,1.5,2.5,4,10,20,50)) 
print(fp_univariate)

# Forest plot multivariate
fp_multi <- ggplot(data=full_model_res, aes(x=var, y=or, ymin=cilo, ymax=ciup)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariate") + ylab("Odd ratio (95% CI)") +
  theme_bw() # use a white background 
fp_multi   <- fp_multi  + scale_y_continuous(trans = log_trans(), breaks = c(0.2,0.5,1,1.5,2.5,4,10,20,50)) 
print(fp_multi )

# Lasso CV lambda
plot(cv_output)




# Environment -------------------------------------------------------------

# R version 3.6.2 (2019-12-12)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19041)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Ireland.1252  LC_CTYPE=English_Ireland.1252    LC_MONETARY=English_Ireland.1252
# [4] LC_NUMERIC=C                     LC_TIME=English_Ireland.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scales_1.1.0   ggplot2_3.3.2  glmnet_4.1     Matrix_1.2-18  mada_0.5.10    mvmeta_1.0.3   ellipse_0.4.2 
# [8] mvtnorm_1.0-11
