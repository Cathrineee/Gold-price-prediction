# import data
data <- read.csv("../data.csv", header=T)[c("Date", "Switzerland.CHF.")]
head(data)
price <- ts(data$Switzerland.CHF., start=c(1979, 1), end=c(2021, 7), frequency=12)
time <- time(price)
boxplot(price)
par(mfrow=c(1,2))
ts.plot(price)
abline(v=2017+1/3, col="red")
acf(price)

decompose.add <- decompose(price, type="additive")
decompose.mult <- decompose(price, type="multiplicative")
plot(decompose.add)
plot(decompose.mult)


# 90% training 10% test set
training.price = window(price, start=c(1979,1), end=c(2017,3), frequency=12)
training.time = window(time, start=c(1979,1), end=c(2017,3), frequency=12)
test.price = window(price, start=c(2017, 4), frequency=12)
test.time = window(time, start=c(2017, 4), frequency=12)
test = data.frame(time=test.time)
training <- data.frame(price=training.price, time=training.time)
df <- data.frame(price=price, time=time)
future.time = time(ts(start=c(2021,7), end=c(2023,12), frequency=12))
future = data.frame(time=future.time)

library(MASS)

## non-constant variance
seg = factor(c(rep(1:4, each=102), rep(5, 103)))
fligner.test(price, seg)

power.alpha <- c(-2,-1.5,-1,-0.5,0,0.5,1.5,2)
bx = boxcox(lm(training.price~1), lambda=power.alpha)
fligner.test(price^bx$x[which.max(bx$y)], seg)
f.p <- c()
for (a in power.alpha) {
  if (a == 0) {
    y.star = log(price)
  } else {
    y.star = price^a
  }
  f.test = fligner.test(y.star, seg)
  f.p = c(f.p, f.test$p.value)
}
f.p
plot(power.alpha, f.p, xlab="Alpha", ylab="P-value")
abline(h=0.1, col="red")
print(power.alpha[which(f.p>0.1)])

power.alpha.opt = bx$x[which.max(bx$y)]
training.price.star = training.price^power.alpha.opt
test.price.star = test.price^power.alpha.opt
price.star = price^power.alpha.opt
fligner.test(training.price.star, seg)
par(mfrow=c(2,1))
ts.plot(training.price.star)
acf(training.price.star)



## regression
degrees <- c(2:10)
MSE <- function(y, yhat) {
  mean((y-yhat)^2)
}
# create fit models with different degrees
poly.mse = c()
for (d in degrees) {
  model = lm(price ~ poly(as.vector(time), d), data=training)
  pred = predict.lm(model, test, interval = "none")
  mse = MSE(test.price, pred)
  poly.mse = c(poly.mse, mse)
}
print(data.frame(degree=degrees, MSE=poly.mse))
poly.opt = degrees[which.min(poly.mse)]
plot(degrees, poly.mse)
print(paste("Degree of the most optimal polynomial is: ", poly.opt))

pd8 = lm(price ~ poly(as.vector(time), 8), data=training)
pd2 = lm(price ~ poly(as.vector(time), 2), data=training)
pd2.pred = predict.lm(pd2, test, interval = "none")
pd8.pred = predict.lm(pd8, test, interval = "none")
ts.plot(price, ylab="Price", main="Polynomial regression with degree 8")
points(training.time, pd8$fitted.values, type="l", col="red")
points(training.time, pd2$fitted.values, type="l", col="blue")
points(test.time, pd8.pred, type="l", col="red")
points(test.time, pd2.pred, type="l", col="blue")
abline(v=2017+1/3, col="green")



library(glmnet)

# possible lambda values
Log.Lambda.Seq = c(seq(-15,15,by=0.5)) 
Lambda.Seq = exp(Log.Lambda.Seq)
set.seed(443)

## RIDGE
lambda.ridge = c()
CV.ridge.error = c()
for (p in degrees) {
  X = poly(as.vector(training.time), degree=p, simple=T)
  CV = cv.glmnet(X, as.vector(training.price), alpha=0, lambda=Lambda.Seq, 
                 standardize = T, family = "gaussian")
  lambda.ridge = c(lambda.ridge, CV$lambda.min)
  CV.ridge.error = c(CV.ridge.error, CV$cvm[CV$index[1,1]])
}

par(mfrow=c(1,2))
plot(degrees, lambda.ridge)
plot(degrees, CV.ridge.error)

degree.ridge.opt = degrees[which.min(CV.ridge.error)]
lambda.ridge.opt = lambda.ridge[which.min(CV.ridge.error)]
sprintf("The best ridge model is at degree %d with lambda %.4f", degree.ridge.opt, lambda.ridge.opt)

## LASSO
lambda.lasso = c()
CV.lasso.error = c()
for (p in degrees) {
  X = poly(as.vector(training.time), degree=p, simple=T)
  CV = cv.glmnet(X, as.vector(training.price), alpha=1, lambda=Lambda.Seq, 
                 standardize = T, family = "gaussian")
  lambda.lasso = c(lambda.lasso, CV$lambda.min)
  CV.lasso.error = c(CV.lasso.error, CV$cvm[CV$index[1,1]])
}

par(mfrow=c(1,2))
plot(degrees, lambda.lasso)
plot(degrees, CV.lasso.error)

degree.lasso.opt = degrees[which.min(CV.lasso.error)]
lambda.lasso.opt = lambda.lasso[which.min(CV.lasso.error)]
sprintf("The best lasso model is at degree %d with lambda %.4f", degree.lasso.opt, lambda.lasso.opt)

## ELASTIC NET
lambda.elastic = c()
CV.elastic.error = c()
for (p in degrees) {
  X = poly(as.vector(training.time), degree=p, simple=T)
  CV = cv.glmnet(X, as.vector(training.price), alpha=0.5, lambda=Lambda.Seq, 
                 standardize = T, family = "gaussian")
  lambda.elastic = c(lambda.elastic, CV$lambda.min)
  CV.elastic.error = c(CV.elastic.error, CV$cvm[CV$index[1,1]])
}

par(mfrow=c(1,2))
plot(degrees, lambda.elastic)
plot(degrees, CV.elastic.error)

degree.elastic.opt = degrees[which.min(CV.elastic.error)]
lambda.elastic.opt = lambda.elastic[which.min(CV.elastic.error)]
sprintf("The best elastic model is at degree %d with lambda %.4f", degree.elastic.opt, lambda.elastic.opt)


ridge.time = poly(as.vector(training.time), degree=degree.ridge.opt, simple=T)
ridge.model = glmnet(ridge.time, as.vector(training.price), alpha=0, lambda=lambda.ridge.opt, 
                     standardize = T, family = "gaussian")
ridge.fit = predict(ridge.model, s=lambda.ridge.opt, 
                    newx=ridge.time,
                    type="response")
lasso.time = poly(as.vector(training.time), degree=degree.lasso.opt, simple=T)
lasso.model = glmnet(lasso.time, as.vector(training.price), alpha=1, lambda=lambda.lasso.opt, 
                     standardize = T, family = "gaussian")
lasso.fit = predict(lasso.model, s=lambda.lasso.opt, 
                     newx=lasso.time, 
                     type="response")
elastic.time = poly(as.vector(training.time), degree=degree.elastic.opt, simple=T)
elastic.model = glmnet(elastic.time, as.vector(training.price), alpha=0.5, lambda=lambda.elastic.opt, 
                       standardize = T, family = "gaussian")
elastic.fit = predict(elastic.model, s=lambda.elastic.opt, 
                       newx=elastic.time, 
                       type="response")
ridge.pred = predict(ridge.model, s=lambda.ridge.opt, 
                     newx=poly(as.vector(test.time), degree=degree.ridge.opt, simple=T), 
                     type="response")
ridge.mse = MSE(test.price, ridge.pred)
lasso.pred = predict(lasso.model, s=lambda.lasso.opt, 
                     newx=poly(as.vector(test.time), degree=degree.lasso.opt, simple=T), 
                     type="response")
lasso.mse = MSE(test.price, lasso.pred)
elastic.pred = predict(elastic.model, s=lambda.elastic.opt, 
                       newx=poly(as.vector(test.time), degree=degree.elastic.opt, simple=T), 
                       type="response")
elastic.mse = MSE(test.price, elastic.pred)

par(mfrow=c(2,2))
ts.plot(price, main="Polynomial Regression", ylim=c(200,2000))
points(training.time, pd2$fitted.values, type="l", col="red")
points(test.time, pd2.pred, type="l", col="blue")
ts.plot(price, main="Ridge Regression", ylim=c(200,2000))
points(training.time, ridge.fit, type="l", col="red")
points(test.time, ridge.pred, type="l", col="blue")
ts.plot(price, main="LASSO Regression", ylim=c(200,2000))
points(training.time, lasso.fit, type="l", col="red")
points(test.time, lasso.pred, type="l", col="blue")
ts.plot(price, main="Elastic Net Regression", ylim=c(200,2000))
points(training.time, elastic.fit, type="l", col="red")
points(test.time, elastic.pred, type="l", col="blue")

all = data.frame(polynomial=MSE(test.price, pd2.pred), ridge=ridge.mse, lasso=lasso.mse, elastic=elastic.mse)
all

poly.model = lm(price ~ poly(as.vector(time), poly.opt), data=df)
summary(poly.model)
poly.pred = predict.lm(poly.model, interval = "none")
poly.MSE = MSE(price, poly.pred)
poly.MSE

## residual diagnostic
res = price - poly.pred
par(mfrow=c(2,2))
plot(poly.pred, res, xlab="Fitted values", ylab="Residuals")
abline(h=0, col="red")
plot(res, ylab="Residuals")
abline(h=0, col="red")
car::qqPlot(res, pch=16, col=adjustcolor("black", 0.7),
            xlab = "Theoretical Quantiles (Normal)",
            ylab = "Sample Quantiles (r.hat)",
            main = "Normal Q-Q Plot")
acf(res, main="Sample ACF of Residuals")

reg.pred = predict(poly.model, future, interval="predict")
ts.plot(price, xlim=c(1979,2024), ylim=c(350,2300),
        main="Polynomial Regression Model Prediction", 
        xlab="Time", ylab="Gold Price")
points(time, poly.pred, col="red",type='l')
polygon(c(future.time, rev(future.time)), c(reg.pred[,3], rev(reg.pred[,2])), col = "grey", border=NA)
points(future.time, reg.pred[,1], col="blue", type="l")


# shapiro-wilk
shapiro.test(res)

# constant variance
seg = c(rep(1:4, each=102), rep(5, 103))
fligner.test(res, seg)

# difference sign test
randtests::difference.sign.test(res)

# run tests
lawstat::runs.test(res)



#### SMOOTHING

# simple exponential smoothing
es<-HoltWinters(training.price, gamma = FALSE, beta=FALSE)
es_pred = predict(es,n.ahead=52)
es_mse <- MSE(test.price, es_pred)

# double exponential smoothing
des<-HoltWinters(training.price, gamma = FALSE)
des_pred = predict(des,n.ahead=52)
des_mse<-MSE(test.price, des_pred)

# hw.additive
hw.additive=HoltWinters(training.price,seasonal = 'additive')
hwa_pred = predict(hw.additive,n.ahead=52)
hwa_mse<-MSE(test.price, hwa_pred)

# hw.multiplicative
hw.multi =HoltWinters(training.price,seasonal = 'multiplicative')
hwm_pred = predict(hw.multi,n.ahead=52)
hwm_mse<-MSE(test.price, hwm_pred)

data.frame(single=es_mse, double=des_mse, additiveHW=hwa_mse, multiplicativeHW=hwm_mse)

par(mfrow=c(2,2))
ts.plot(price, main="Single Exponential Smoothing")
points(training.time[2:length(training.time)], es$fitted[,1], type="l", col="red")
points(test.time, es_pred, type="l", col="blue")
ts.plot(price, main="Double Exponential Smoothing")
points(training.time[3:length(training.time)], des$fitted[,1], type="l", col="red")
points(test.time, des_pred, type="l", col="blue")
ts.plot(price, main="Additive Holt-Winters")
points(training.time[13:length(training.time)], hw.additive$fitted[,1], type="l", col="red")
points(test.time, hwa_pred, type="l", col="blue")
ts.plot(price, main="Multiplicative Holt-Winters")
points(training.time[13:length(training.time)], hw.multi$fitted[,1], type="l", col="red")
points(test.time, hwm_pred, type="l", col="blue")

# double exponential smoothing - entire data-set
smooth <- HoltWinters(price, gamma = FALSE)
smooth.fit = smooth$fitted[,1]
smooth.mse = MSE(price[3:length(price)], smooth.fit)
smooth_pred = predict(smooth, n.ahead=30, prediction.interval = T)

## residual diagnostic
res = price - smooth.fit
par(mfrow=c(2,2))
plot(smooth.fit, res, xlab="Fitted values", ylab="Residuals")
abline(h=0, col="red")
plot(res, ylab="Residuals")
abline(h=0, col="red")
car::qqPlot(res, pch=16, col=adjustcolor("black", 0.7),
            xlab = "Theoretical Quantiles (Normal)",
            ylab = "Sample Quantiles (r.hat)",
            main = "Normal Q-Q Plot")
acf(res, main="Sample ACF of Residuals")

# shapiro-wilk
shapiro.test(res)

# constant variance
seg = c(rep(1:4, each=102), rep(5, 103))
fligner.test(res, seg)

# difference sign test
randtests::difference.sign.test(res)

# run tests
lawstat::runs.test(res)

ts.plot(price, xlim=c(1979,2024), ylim=c(350,2300),
        main="Double Exponential Smoothing Method Prediction", 
        xlab="Time", ylab="Gold Price")
points(time[3:length(time)], smooth.fit, col="red",type='l')
polygon(c(future.time, rev(future.time)), c(smooth_pred[,3], rev(smooth_pred[,2])), col = "grey", border=NA)
points(future.time, smooth_pred[,1], col="blue", type="l")





## BOX-JENKINS

# one regular differencing
diff1 = diff(training.price)
par(mfrow=c(1,3))
plot(diff1, ylab="Differenced Price", main="First-order Regular Differencing")
acf(diff1)
pacf(diff1)
# second regular differencing...?
diff2 = diff(diff1)
par(mfrow=c(1,3))
plot(diff2, ylab="Differenced Price", main="Second-order Regular Differencing")
acf(diff2)
pacf(diff2)

# fit models
library(astsa)
model1 = sarima(training.price, 0,1,1)
model2 = sarima(training.price, 1,1,0)
model3 = sarima(training.price, 0,1,5)
model4 = sarima(training.price, 0,1,6)
model5 = sarima(training.price, 6,1,0)
model6 = sarima(training.price, 1,1,1)
model7 = sarima(training.price, 1,1,2)
model8 = sarima(training.price, 2,1,1)
model9 = sarima(training.price, 2,1,2)

pred_model1 = sarima.for(training.price, p=0, d=1, q=1, n.ahead = 52)
pred_model2 = sarima.for(training.price, p=1, d=1, q=0, n.ahead = 52)
pred_model3 = sarima.for(training.price, 0,1,5, n.ahead = 52)
pred_model4 = sarima.for(training.price, 0,1,6, n.ahead = 52)
pred_model5 = sarima.for(training.price, 6,1,0, n.ahead = 52)
pred_model6 = sarima.for(training.price, p=1, d=1, q=1, n.ahead = 52)
pred_model7 = sarima.for(training.price, 1,1,2, n.ahead = 52)
pred_model8 = sarima.for(training.price, 2,1,1, n.ahead = 52)
pred_model9 = sarima.for(training.price, 2,1,2, n.ahead = 52)
ds1 = data.frame(
  Model = c("ARIMA(0,1,1)", "ARIMA(1,1,0)", "ARIMA(0,1,5)",
            "ARIMA(0,1,6)", "ARIMA(6,1,0)", "ARIMA(1,1,1)",
            "ARIMA(1,1,2)", "ARIMA(2,1,1)", "ARIMA(2,1,2)"),
  MSE = c(MSE(test.price, pred_model1$pred),MSE(test.price, pred_model2$pred),MSE(test.price, pred_model3$pred),
          MSE(test.price, pred_model4$pred),MSE(test.price, pred_model5$pred),MSE(test.price, pred_model6$pred),
          MSE(test.price, pred_model7$pred),MSE(test.price, pred_model8$pred),MSE(test.price, pred_model9$pred)),
  AIC = c(model1$AIC, model2$AIC, model3$AIC,
          model4$AIC, model5$AIC, model6$AIC,
          model7$AIC, model8$AIC, model9$AIC),
  AICc = c(model1$AICc, model2$AICc, model3$AICc,
           model4$AICc, model5$AICc, model6$AICc,
           model7$AICc, model8$AICc, model9$AICc),
  BIC = c(model1$BIC,model2$BIC,model3$BIC,
          model4$BIC,model5$BIC,model6$BIC,
          model7$BIC,model8$BIC,model9$BIC)
)
ds1

# residual diagnostic of ARIMA(6,1,0)
fit_train5 = training.price - model5$fit$residuals
fit_train1 = training.price - model1$fit$residuals
fit_train2 = training.price - model2$fit$residuals
par(mfrow=c(1,3))
ts.plot(price, main="ARIMA(0,1,1) fit on Training")
points(training.time, fit_train1, type="l", col="red")
points(test.time, pred_model1$pred, type="l", col="blue")
ts.plot(price, main="ARIMA(1,1,0) fit on Training")
points(training.time, fit_train2, type="l", col="red")
points(test.time, pred_model2$pred, type="l", col="blue")
ts.plot(price, main="ARIMA(6,1,0) fit on Training")
points(training.time, fit_train5, type="l", col="red")
points(test.time, pred_model5$pred, type="l", col="blue")

par(mfrow=c(2,2))
plot(fit_train, model5$fit$residuals, xlab="Fitted values", ylab="Residuals")
abline(h=0, col="red")
plot(model5$fit$residuals, ylab="Residuals")
abline(h=0, col="red")
car::qqPlot(model5$fit$residuals, pch=16, col=adjustcolor("black", 0.7),
            xlab = "Theoretical Quantiles (Normal)",
            ylab = "Sample Quantiles (r.hat)",
            main = "Normal Q-Q Plot")
acf(model6$fit$residuals, main="Sample ACF of Residuals")

#assumption test
shapiro.test(model5$fit$residuals)
randtests::difference.sign.test(model5$fit$residuals)

# fit to all data
final_arima = sarima(price, 6,1,0)
fit_val = price - final_arima$fit$residuals

# MSE, AIC, BIC
ds2 = data.frame(
  Model = c("ARIMA(6,1,0)"),
  MSE = c(MSE(price, fit_val)),
  AIC = c(final_arima$AIC),
  AICc = c(final_arima$AICc),
  BIC = c(final_arima$BIC)
)
ds2

pred = sarima.for(price,6,1,0, n.ahead = 29, main="ARIMA(6,1,0) Prediction")
lower <- pred$pred-1.96*pred$se
upper <- pred$pred+1.96*pred$se
x = c(time(upper) , rev(time(upper)))
y = c(upper , rev(lower))
plot(price, xlim = c(1979,2024), ylim = c(370, 2100))
polygon(x, y, col = "grey" , border =NA)
lines(pred$pred, col = "red")


## final comparison
data.frame(model=c("Polynomial", "Smoothing", "Box-Jenkins"),
           MSE.fit = c(poly.MSE, smooth.mse, MSE(price, fit_val)),
           MSE.pred = c(MSE(test.price, pd2.pred), des_mse, MSE(test.price, pred_model5$pred)))



