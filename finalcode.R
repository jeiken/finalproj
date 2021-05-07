library(glmnet)

## Data wrangling and normalization 

load("DOC_baseFlow_weightedAvg_and_predictors.RData")

raw_data <- DOC_baseFlow_weightedAvg_and_predictors
y_response <- raw_data[,"meanDOC"]


xFactors <- model.matrix(y_response ~ as.factor(raw_data[,"aspect_predominant"]))[, -1]

colNames_xFactors <- c("aspect_predominant_1")
for (index in 2:8)
{
  colNames_xFactors <- c(colNames_xFactors,paste("aspect_predominant_",as.character(index),sep=""))
}

colnames(xFactors) <- colNames_xFactors

X_predictors_unStandardized <- as.matrix( cbind( DOC_baseFlow_weightedAvg_and_predictors[,c(2,6,9:30,32)], xFactors) )
dataFrame_X_predictors <- as.data.frame(X_predictors_unStandardized)

X_predictors_means <- colMeans(X_predictors_unStandardized)
X_predictors_sd <- apply(X_predictors_unStandardized, MARGIN=2, 'sd')

time_index <- 			1
fire_indices <- 		2
MPB_indecies <- 		3:5
landCover_indices <- 	6:8
temp_indices <- 		9:14
precipSnow_indices <- 	15:20
soil_indices <- 		21:23
elevation_index <- 		24
wasteWater_index <- 	25
aspec_indices <- 		26:33
X_indices_subtractMean <- c(time_index, temp_indices, precipSnow_indices, soil_indices, elevation_index, wasteWater_index)
X_indices_scaleBySD <- c(fire_indices)

X_predictors_standardized <- X_predictors_unStandardized
X_predictors_standardized[,X_indices_subtractMean] <- sweep(X_predictors_unStandardized[,X_indices_subtractMean],MARGIN=2,X_predictors_means[X_indices_subtractMean],'-')

colnames(X_predictors_standardized)
data <- dataFrame_X_standardized <- as.data.frame(X_predictors_standardized)


### Model Creation 

set.seed(222)

n <- nrow(data)
train_size <- 939  # 80% of the data (20% left for testing)

backlist_r2 <- c()
forlist_r2 <- c()
lassolist_r2 <- c()
backlist_RMSE <- c()
forlist_RMSE <- c()
lassolist_RMSE <- c()
backlist_MAE <- c()
forlist_MAE <- c()
lassolist_MAE <- c()


for( i in 1:10){
  p <- sample(1:n)
  x_train <- data[p[1:train_size], ]
  x_test <- data[p[(train_size + 1):n], ]
  y_train <- y_response[p[1:train_size]]
  y_test <- y_response[p[(train_size + 1):n]]
  train <- cbind(y_train, x_train)
  test <- cbind(y_test, x_test)
  mean <- mean(train$y_train)
  
  ## OG model
  model <- lm(formula = y_train ~ ., data = train)
  summary(model)
  
  
  ## Backwards 
  BIC <- step(model, direction = "backward", k = log(nrow(train)), trace = 0)
  summary(BIC) # 11/31 variables selected, R^2 = 0.6154
  
  predictions_BIC <- predict(BIC, test)
  RMSE_BIC <- sqrt(mean(((predictions_BIC - y_test)^2)))
  MAE_BIC <- mean(abs(predictions_BIC - y_test))
  
  BIC_residuals <- predictions_BIC - y_test
  BIC_SSR <- sum(BIC_residuals^2) 
  
  SST <- sum((y_test - mean)^2) # 7336.643
  BIC_rSquared <- 1 - BIC_SSR/SST

  backlist_r2 <- c(backlist_r2, BIC_rSquared)
  backlist_RMSE <- c(backlist_RMSE, RMSE_BIC)
  backlist_MAE <- c(backlist_MAE, MAE_BIC)
  
  
  ## Forwards
  small_model <- lm(y_train ~ 1, data = train)
  BICF <- step(small_model, direction = "forward", scope = formula(model), k = log(nrow(train)), trace = 0)
  summary(BICF) # 9 variables selected, R^2 = 0.6119
  
  predictions_BICF <- predict(BICF, test)
  RMSE_BICF <- sqrt(mean(((predictions_BICF - y_test)^2)))
  MAE_BICF <- mean(abs(predictions_BICF - y_test))
  
  BICF_residuals <- predictions_BICF - y_test
  BICF_SSR <- sum(BICF_residuals^2) 
  
  SST <- sum((y_test - mean)^2) # 7336.643
  BICF_rSquared <- 1 - BICF_SSR/SST
  
  forlist_r2 <- c(forlist_r2, BICF_rSquared)
  forlist_RMSE <- c(forlist_RMSE, RMSE_BICF)
  forlist_MAE <- c(forlist_MAE, MAE_BICF)
  
  ## Lasso 
  leave_out <- which((colnames(train) == "y_train"))
  x <- as.matrix(train[,-leave_out])
  y <- train[, "y_train"]
  lasso <- cv.glmnet(x, y)
  coef(lasso) # selected 5 variables,
  
  leave_out2 <- which((colnames(test) == "y_test"))
  x2 <- as.matrix(test[,-leave_out2])
  y2 <- test[, "y_test"]
  
  lasso_predictions <- predict(lasso, newx = x2)
  lasso_residuals <- lasso_predictions - y_test
  lasso_SSR <- sum(lasso_residuals^2) # 3368.946
  SST <- sum((y_test - mean)^2) # 7336.643
  lasso_rSquared <- 1 - lasso_SSR/SST #unadjusted 0.5408055
  
  RMSE_lasso <- sqrt(mean(((lasso_predictions - y_test)^2)))
  MAE_lasso <- mean(abs(lasso_predictions - y_test))
  
  lassolist_r2 <- c(lassolist_r2, lasso_rSquared)
  lassolist_RMSE <- c(lassolist_RMSE, RMSE_lasso)
  lassolist_MAE <- c(lassolist_MAE, MAE_lasso)
  
  print(i)
  
}

mean(backlist_r2)
mean(forlist_r2)
mean(lassolist_r2)
mean(backlist_RMSE)
mean(forlist_RMSE)
mean(lassolist_RMSE)
mean(backlist_MAE)
mean(forlist_MAE)
mean(lassolist_MAE)


#notes: 
# use lasso to select variables, then fit with lm
# how many variables each method selects in common with the others

## Model Creation pt 2

coef(lasso)

new_model <- lm(y_train ~ maxTemp_baseFlow + maxTemp_yearAvg + soil_Afraction + elevation_mean + wasteWaterPointSources_count + aspect_predominant_3, data = train)
summary(new_model)


summary(BIC)
summary(BICF)


newest_model <- lm(y_train ~ elevation_mean + wasteWaterPointSources_count + aspect_predominant_3, data = train)
summary(newest_model)
