#Basic examples
library(seasonality)

#this can also be entered as a ts_boxable object
dates_weekly <- as.Date(initial_claims$DATE)
data_weekly  <- initial_claims$ICNSA

SA <- auto_SA(data = data_weekly, dates = dates_weekly)

ts.plot(cbind(data_weekly, SA$data_sa), col = c("red", "blue")) #could use some fine tuining but pretty good!

dt <- initial_claims

monthly <- agg_to_monthly(dt)

data_monthly <- monthly$ICNSA
dates_monthly <- as.Date(monthly$DATE)

SA <- auto_SA(data = data_monthly, dates = dates_monthly)

# (1) Do factors add up?
zero <- data_monthly - SA$data_sa - SA$adj_factor
max(zero) #yes

#(2) test against other packages (other stuff works for monthly data)
##### test against the forecast package ######################

library(forecast)
data_ts <- ts(data_monthly, start = c(1967,1), frequency = 12)

decomposed <- decompose(data_ts)
sa <- forecast::seasadj(decomposed)

#Compare results

ts.plot(cbind(sa, SA$data_sa), col = c("red", "blue")) #no contest

##### test against the seasonal package ######################

library(seasonal)
sa <- seas(data_ts, x11 = "", transform.function = "none") #, x11.appendfcst = "yes"
sa_factor <- series(sa, "x11.adjustfac") #test the factors
zero <- data_ts - final(sa) - sa_factor
max(zero) #check
ts.plot(cbind(final(sa), SA$data_sa), col = c("red", "blue")) #very similar... can't say who's best by just looking
#Which series offers better predictive power?
first <- auto.arima(SA$data_sa[-(1:13)])
second <- auto.arima(final(sa)[-(1:13)])
print(first)
print(second)

#seasonal does better, but at least we're close!



