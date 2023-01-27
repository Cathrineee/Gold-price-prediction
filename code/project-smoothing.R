
price <- ts(data$Switzerland.CHF., start=c(1979, 1), end=c(2021, 7), frequency=12)

plot(price)

# Flingner 

seg<-c(rep(1:5,each=102),1)
seg
alpha <- c(-1,-0.75,-0.5,0,0.5,2,3)
alpha.opt <- 0
f.p <- 0
for (a in alpha) {
  if (a == 0) {
    y.star = log(price)
  } else {
    y.star = price^a
  }
  f.test = fligner.test(y.star, seg)
  if (f.test$p.value > f.p) {
    alpha.opt = a
    f.p = f.test$p.value
  }
}
print(paste("Alpha of the optimal transformation: ", alpha.opt))

fligner.test(price^(-1), seg)
plot(price^(-1))

# just to see how the data works
a<-decompose(price)
plot(a)

# split for training and test
price_train<-window(price, 1979, 2012.999)
price_train
price_test<-window(price, 2013, c(2021,7))
price_test

# another way
data <- read.csv("Data_Group2.csv", header=T)[c("Date", "Switzerland.CHF.")]
price <- ts(data$Switzerland.CHF., start=c(1979, 1), end=c(2021, 7), frequency=12)
time <- time(price)


test.idx = sample(1:length(price), 100)
test.price <- price[test.idx]
price_test<- data.frame(time=time[test.idx])

training.price <- price[-test.idx]
training.time <- time[-test.idx]
price_train <- data.frame(price=training.price, time=training.time)

# simple exponential smoothing
es<-HoltWinters(price_train, gamma = FALSE, beta=FALSE)
es
es_pred = predict(es,n.ahead=12)
es_pred
es_mse<-mean((price_test-es_pred)^2)
es_mse

plot(price)
plot(es)
abline(es, col="red")

# double exponential smoothing
des<-HoltWinters(price_train, gamma = FALSE)
des
des_pred = predict(des,n.ahead=12)
des_pred
des_mse<-mean((price_test-des_pred)^2)
des_mse

# hw.additive
hw.additive=HoltWinters(price_train,seasonal = 'additive')
hwa_pred = predict(hw.additive,n.ahead=12)
hwa_pred
hwa_mse<-mean((price_test-hwa_pred)^2)
hwa_mse

# hw.multiplicative
hw.multi =HoltWinters(price_train,seasonal = 'multiplicative')
hwm_pred = predict(hw.multi,n.ahead=12)
hwm_pred
hwm_mse<-mean((price_test-hwm_pred)^2)
hwm_mse

c(es_mse, des_mse, hwa_mse, hwm_mse)      
                  
# best model: simple exponential

best<-HoltWinters(price, gamma = FALSE, beta=FALSE)
res<- price-best$fitted[,1]
plot(res)
acf(res)
shapiro.test(res)

plot(seq(1,512), res)


# test for regression

training.price <- window(price, start = c(1979,1), end = c(2012,12))
training.time <- as.vector(time(training.price))
training.time
test.price <- window(price, start = c(2013,1), end = c(2021,7))
new <- data.frame(time=seq(2013, 2021, length=103))
new

d=4
model = lm(price_train ~ poly(training.time, d))
pred = predict.lm(model, new, interval = "none")
pred





