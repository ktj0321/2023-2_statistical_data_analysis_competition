## 라이브러리 추가
library(lmtest)
library(TSA)
library(fpp)
library(tseries)

## 자료 읽기
data <- read.csv(file.choose(), header = TRUE)
attach(data)


################################################################################
####################### Arima 오차항을 갖는 회귀모형  ##########################
## 암 자료 적합
plot(cancer, type = "l") # 증가하는 추세가 보이므로, 차분
cancerd <- diff(cancer)
plot(cancerd, type = "l") # 정상적 시계열로 보임.
acf(cancerd)
pacf(cancerd)
# ARIMA 모형이 식별되지 않으므로, 회귀모형 적합을 생각해볼 수 있음.

# 잔차가 독립인지 검정 (Durbin-Watson 검정)
dwtest(cancer~yts) # 귀무가설 기각. 잔차들의 자기상관이 존재함.

# 잔차들이 Arima 모형을 따르는지 알아보기 위해,
canfit <- lm(cancer ~ yts) # 단순 선형회귀모형 적합.
summary(canfit)
rs_canfit <- residuals(canfit) # 해당 모형의 잔차
plot(rs_canfit) # 잔차의 시퀀스 차트. 패턴을 보이는 것처럼 보임.
abline(h=0, col='blue')
acf(rs_canfit)
pacf(rs_canfit) # 잔차가 AR(1)으로 식별됨.

# Arima 오차항을 갖는 회귀모형, ARIMA(1, 0, 0)
canfit1 <- Arima(cancer, xreg=yts, order=c(1, 0, 0))
canfit1


summary(canfit1)

# 진단
rs_canfit1 <- residuals(canfit1)
Box.test(rs_canfit1, fitdf=1, lag=20, type="Ljung")

plot(rs_canfit1, type = "p")
abline(h=0, col='blue')
acf(rs_canfit1)
pacf(rs_canfit1)

# 암에 의한 사망자 수 향후 5년간 예측 (index 39는 2023년을 의미)
fcast_cancer <- forecast(canfit1, xreg=39:44, h=5)
fcast_cancer

# 그리기
plot(cancer, type = "l", lwd = 3)
lines(fitted(canfit1), type = "l", lwd = 3, col = "red")
legend("topleft", legend = c("관측값(원래 시계열 자료)", "모형의 적합값"), # 범례 추가
       col = c("black", "red"), lty = c(1, 1))




## 신경질환 자료 적합
plot(neu, type = "l") # 급격하게 증가하는 경향 보이므로, log 취함.
logneu <- log(neu)
plot(logneu, type = 'l') # 증가하는 추세를 보이므로, 차분
logneud <- diff(logneu)
plot(logneud, type = "l") # 정상적 시계열로 보임.
acf(logneud)
pacf(logneud)
# ARIMA 모형이 식별되지 않으므로, 회귀모형 적합을 생각해볼 수 있음.

# 잔차가 독립인지 검정 (Durbin-Watson 검정)
dwtest(logneu~yts) # 귀무가설 기각. 잔차들의 자기상관이 존재함.

# 잔차들이 Arima 모형을 따르는지 알아보기 위해,
neufit <- lm(logneu ~ yts) # 단순회귀모형 적합해봄.
rs_neufit <- residuals(neufit) # 해당 모형의 잔차
plot(rs_neufit) # 잔차의 시퀀스 차트. 패턴을 보이는 것처럼 보임.
abline(h=0, col='blue')
acf(rs_neufit)
pacf(rs_neufit) # 잔차가 AR(1)으로 식별됨.

# Arima 오차항을 갖는 회귀모형, ARIMA(1, 0, 0)
neufit1 <- Arima(logneu, xreg=yts, order=c(1, 0, 0))
neufit1

# 진단
rs_neufit1 <- residuals(neufit1)
Box.test(rs_neufit1, fitdf=1, lag=20, type="Ljung")

plot(rs_neufit1)
abline(h=0, col='blue')
acf(rs_neufit1)
pacf(rs_neufit1)

# 예측
fcast_neu <- forecast(neufit1, xreg=39:44, h=5)
fcast_neu

# 그리기
plot(neu, type = "l", lwd = 3)
lines(exp(fitted(neufit1)), type = "l", lwd = 3, col = "red")
legend("topleft", legend = c("관측값(원래 시계열 자료)", "모형의 적합값"), # 범례 추가
       col = c("black", "red"), lty = c(1, 1))



## 심혈관질환 자료 적합
plot(cyl, type = "l") # 감소 추세가 보이므로, 차분
cyld <- diff(cyl)
plot(cyld, type = 'l') # 정상적 시계열로 보임.
acf(cyld)
pacf(cyld)
# ARIMA 모형이 식별되지 않으므로, 회귀모형 적합을 생각해볼 수 있음.

# 잔차가 독립인지 검정 (Durbin-Watson 검정)
dwtest(cyl~yts) # 귀무가설 기각. 잔차들의 자기상관이 존재함.

# 잔차들이 Arima 모형을 따르는지 알아보기 위해,
cylfit <- lm(cyl ~ yts) # 단순회귀모형 적합해봄.
rs_cylfit <- residuals(cylfit) # 해당 모형의 잔차
plot(rs_cylfit) # 잔차의 시퀀스 차트. 패턴을 보이는 것처럼 보임.
abline(h=0, col='blue')
acf(rs_cylfit) # 잔차가 MA(2)로 식별됨.
pacf(rs_cylfit) # 잔차가 AR(1)으로 식별됨.

# Arima 오차항을 갖는 회귀모형, ARIMA(1, 0, 0)
cylfit1 <- Arima(cyl, xreg=yts, order=c(1, 0, 0))
cylfit1

# 진단
rs_cylfit1 <- residuals(cylfit1)
Box.test(rs_cylfit1, fitdf=1, lag=20, type="Ljung")

plot(rs_cylfit1)
abline(h=0, col='blue')
acf(rs_cylfit1)
pacf(rs_cylfit1)

# 예측
fcast_cyl <- forecast(cylfit1, xreg=39:44, h=5)
fcast_cyl

# 그리기
plot(cyl, type = "l", lwd = 3)
lines(fitted(cylfit1), type = "l", lwd = 3, col = "red")
legend("topright", legend = c("관측값(원래 시계열 자료)", "모형의 적합값"), # 범례 추가
       col = c("black", "red"), lty = c(1, 1))




## 자살로 인한 사망 적합
plot(sui, type = "l") # 증가 추세를 보이므로, 차분
suid <- diff(sui)
plot(suid, type = "l") # 정상적 시계열로 보임.
acf(suid)
pacf(suid)
# ARIMA 모형이 식별되지 않으므로, 회귀모형 적합을 생각해볼 수 있음.

# 잔차가 독립인지 검정 (Durbin-Watson 검정)
dwtest(sui~yts) # 귀무가설 기각. 잔차들의 자기상관이 존재함.

# 잔차들이 Arima 모형을 따르는지 알아보기 위해,
suifit <- lm(sui ~ yts) # 단순회귀모형 적합해봄.
rs_suifit <- residuals(suifit) # 해당 모형의 잔차
plot(rs_suifit) # 잔차의 시퀀스 차트. 패턴을 보이는 것처럼 보임.
abline(h=0, col='blue')
acf(rs_suifit)
pacf(rs_suifit) # 잔차가 AR(1)으로 식별됨.

# Arima 오차항을 갖는 회귀모형, ARIMA(1, 0, 0)
suifit1 <- Arima(sui, xreg=yts, order=c(1, 0, 0))
suifit1

# 진단
rs_suifit1 <- residuals(suifit1)
Box.test(rs_suifit1, fitdf=1, lag=20, type="Ljung") # 잔차들이 독립. 문제 X

plot(rs_suifit1)
abline(h=0, col='blue')
acf(rs_suifit1)
pacf(rs_suifit1)

# 예측
fcast_sui <- forecast(suifit1, xreg=39:44, h=5)
fcast_sui

# 그리기
plot(sui, type = "l", lwd = 3)
lines(fitted(suifit1), type = "l", lwd = 3, col = "red")
legend("topleft", legend = c("관측값(원래 시계열 자료)", "모형의 적합값"), # 범례 추가
       col = c("black", "red"), lty = c(1, 1))



################################################################################
##################################  ARIMA  #####################################
## 호흡기질환 자료 적합
plot(bre, type = "l") # 증가 추세가 보이므로, 차분
bred <- diff(bre)
plot(bred, type = "l") # 여전히 증가 추세가 보이므로, 한 번 더 차분.
bredd <- diff(bred)
plot(bredd, type = "l") # 정상적 시계열로 보임.
acf(bredd) # MA(1)으로 식별됨.
pacf(bredd) # AR(1)으로 식별됨.

# 두 개의 ARIMA 모형을 적합하고, 둘 중 더 적절한 모형 선정.

# ARIMA(0, 2, 1) 적합
brefit2 <- arima(bre, order = c(0, 2, 1))
brefit2

AIC(brefit2) # AIC 값이 189.92.

# ARIMA(1, 2, 0) 적합
brefit1 <- arima(bre, order = c(1, 2, 0))
AIC(brefit1) # AIC 값이 197.5896

# ARIMA(0, 2, 1)을 적합하는 것이 더욱 적절함.
# 진단
rs_brefit <- residuals(brefit2)
Box.test(rs_brefit, lag = 20, type = "Ljung-Box", fitdf = 1)

# 호흡계 질환에 의한 사망자 수 향후 5년 예측
predict_bre <- predict(brefit2, n.ahead = 5)
predict_bre

# 그리기
plot(bre, type = "l", lwd = 3)
lines(as.vector(fitted(brefit2)), type = "l", lwd = 3, col = "red")
legend("topleft", legend = c("관측값(원래 시계열 자료)", "모형의 적합값"), # 범례 추가
       col = c("black", "red"), lty = c(1, 1))



## 소화기질환 자료 적합.
plot(dig, type = "l") # 감소 추세가 보이므로, 차분
digd <- diff(dig)
plot(digd, type = "l") # 증가 추세가 보이므로, 한 번 더 차분
digdd <- diff(digd)
plot(digdd, type = "l") # 정상적 시계열로 보임.
acf(digdd) # MA(1)으로 식별됨.
pacf(digdd) # AR(1)으로 식별됨.

# 두 개의 ARIMA 모형을 적합하여, 더 적합한 모형 결정.
# ARIMA(1, 2, 0) 적합
digfit1 <- arima(dig, order = c(1, 2, 0))
AIC(digfit1) # AIC 값이 124.8669

# ARIMA(0, 2, 1) 적합
digfit2 <- arima(dig, order = c(0, 2, 1))
digfit2
AIC(digfit2) # AIC 값이 123.1049

# ARIMA(0, 2, 1)을 적합하는 것이 더욱 적절함.
# 진단
rs_digfit <- residuals(digfit2)
Box.test(rs_digfit, lag = 20, type = "Ljung-Box", fitdf = 1)

# 예측
predict(digfit2, n.ahead = 5)

# 그리기
plot(dig, type = "l", lwd = 3)
lines(fitted(digfit2), type = "l", lwd = 3, col = "red")
legend("topright", legend = c("관측값(원래 시계열 자료)", "모형의 적합값"), # 범례 추가
       col = c("black", "red"), lty = c(1, 1))



## 운수사고 사망 적합
plot(car_acc, type = "l") # 감소 추세가 보이므로, 차분
car_accd <- diff(car_acc)
plot(car_accd, type = "l") # 정상적 시계열로 보임.
pp.test(car_accd) # 정상적 시계열임을 확인
acf(car_accd) # AR(1)으로 식별됨
pacf(car_accd) # MA(1)으로 식별됨

# 두 개의 ARIMA 모형을 적합하여, 더 적합한 모형 결정.
# ARIMA(1, 1, 0) 적합
car_accfit1 <- arima(car_acc, order = c(1, 1, 0))
AIC(car_accfit1) # AIC 값이 177.835
car_accfit1

# ARIMA(0, 1, 1) 적합
car_accfit2 <- arima(car_acc, order = c(0, 1, 1))
AIC(car_accfit2) # AIC 값이 179.6614

# ARIMA(1, 1, 0)을 적합하는 것이 더욱 적절함.
# 진단
rs_car_accfit <- residuals(car_accfit1)
Box.test(rs_car_accfit, lag = 20, type = "Ljung-Box", fitdf = 1)

# 예측
predict(car_accfit1, n.ahead = 5)

# 그리기
plot(car_acc, type = "l", lwd = 3)
lines(fitted(car_accfit1), type = "l", lwd = 3, col = "red")
legend("topright", legend = c("관측값(원래 시계열 자료)", "모형의 적합값"), # 범례 추가
       col = c("black", "red"), lty = c(1, 1))

