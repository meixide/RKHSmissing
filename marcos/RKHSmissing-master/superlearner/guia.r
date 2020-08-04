library("SuperLearner")
library("ranger")
library("glmnet")
library("randomForest")
library("ggplot2")
library("mgcv")
library("gam")
library("kernlab")


data(Boston, package = "MASS")
# https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html
# https://watermark.silverchair.com/kww165.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAo0wggKJBgkqhkiG9w0BBwagggJ6MIICdgIBADCCAm8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMrGxkLM5ch6TxTBV-AgEQgIICQOxxNOf6htZLTwfL80K4YAShBEE5rjSFXv4iVQyBw0g0S8TBl1decALxzJXPaEG0MY1lQeMjVBodNx5liME3k3eAgAo1f5UEr7T-JXHVroK2hf-PZZzZZatW5t_hxwwddAclzwuV6X39Bd7WZDcLAIFdPPni4_ASm3S_ABSdLMVC-HzGA90sHMoJjG665nvTS9eEaD8sfwRm6MWLJc-NDPwTYwuVbGYiHYcvoW4RupF7y6jtbhz5QtsxWHyN5C_EfAQCvMj2eD_0zF73qULnnu5qNQkTXPfs_p2uMigAUJgybT3RsgzR3SRs3B6S60njrCfzyorDVA80ZPwyt-yvFtGUzsbChDklvtVW6IC_X9hyxG_GwJO89XUBzty1LhKnSad1ZyjAC9tsa5QIV1qtLM6jPEs5LILSJWhb27yTgVh_vil6lBSTP6rvU-XekukdHFHNBGUdux4hgsk69PhMvLISCV89DU5pe_ou2CqCQgu1DpuW4yOgaDXXPSm7MY7GL6qHT0cKbGjNtnqrIJH32e-q7LiO7VnnHv0Bg-jgQ11DQmeNUYDzxUZAwPQ0jYG9N2nbGw2wCLzu8As2FiHPAzxk_qOsJx02Voxhhw1ArBU9w-5mvAJH11awZvOEdRQYeH3yHwKoZgA3BNjFuoQwZfNHluPv8Ao6YXgLlu534k_BKkJp7xheCApX5RgvUNj0ByVAp3RCC0XJ4RMH32hMxQKjxJhHWBT4QIMxmB0IY1eOSfDvTyvwoXNxkWy6acAq7w
outcome = Boston$medv

# Create a dataframe to contain our explanatory variables.
data = subset(Boston, select = -medv)

# Check structure of our dataframe.
str(data)



train_obs = sample(nrow(data), 150)

# X is our training sample.
x_train = data[train_obs, ]

# Create a holdout set for evaluating model performance.
# Note: cross-validation is even better than a single holdout sample.
x_holdout = data[-train_obs, ]

# Create a binary outcome variable: towns in which median home value is > 22,000.
outcome_bin = as.numeric(outcome > 22)

y_train = outcome_bin[train_obs]
y_holdout = outcome_bin[-train_obs]





# This will take about 2x as long as the previous SuperLearner.
cv_sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
                        # For a real analysis we would use V = 10.
                        V = 10,SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"))




cv_sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
                        # For a real analysis we would use V = 10.
                        V = 10,SL.library = c("All"))






cv_sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
                        # For a real analysis we would use V = 10.
                        V = 10,SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.gam", "SL.ksvm", "SL.knn"))







plot(cv_sl) + theme_bw()

summary(cv_sl)

cv_sl

sl= cv_sl







sl = SuperLearner(Y = y_train, X = x_train, family = binomial(),
                  # For a real analysis we would use V = 10.,
                  SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.gam", "SL.ksvm", "SL.knn"))
class(x_holdout)


predict(sl)




pred = predict(sl, x_holdout, onlySL = TRUE)