#Basic examples
library(seasonality)
#
# library(data.table)
# load("~/seasonality/data/icnsa.RData")
# load("~/seasonality/data/holiday2.RData")


# library(devtools)
# load_all("~/seasonality")

# initial_claims_SA <- initial_claims_SA[-nrow(initial_claims_SA),]
# save(initial_claims_SA, file = "~/seasonality/data/icsa.RData")

#this can also be entered as a ts_boxable object
dates <- as.Date(initial_claims$DATE)
data  <- initial_claims$ICNSA

start_time <- Sys.time()
SA <- auto_SA(data = data, dates = dates)
end_time <- Sys.time()
end_time - start_time

#ts.plot(cbind(data_weekly, SA$data_sa), col = c("red", "blue")) #could use some fine tuining but pretty good!
plot(dates, data, type = "l", lwd = 3, col = "lightskyblue1",
     xlab = "Dates", ylab = "Initial Claims")
lines(dates, SA$data_sa, col = "black")

#ts.plot(cbind(data_weekly, SA$data_sa), col = c("red", "blue")) #could use some fine tuining but pretty good!
plot(dates, SA$data_sa, type = "l", lwd = 3, col = "lightskyblue1",
     xlab = "Dates", ylab = "Initial Claims")
lines(dates, initial_claims_SA$ICSA, col = "black")







data <- IP$IPNSA
dates <- IP$DATE

SA <- auto_SA(data = data, dates = dates)

y <- ts(data, start = c(1970,1), frequency = 12)
library(seasonal)

sa <- seas(y)
plot(dates[-1], diff(log(SA$data_sa)), type = "l", lwd = 1, col = "blue",
     xlab = "Dates", ylab = "IP")
lines(dates[-1], diff(log(final(sa))), lwd = 1, col = "red")


#test agains auto.arima
library(forecast)

data <- IP$IPNSA
dates <- IP$DATE

trend <- loess(data ~ seq(length(data)), span = .1, na.action = na.exclude) #refine parameter selection
trend <- predict(trend)
#ts.plot(cbind(trend, data), col = c("red", "blue"))
y <- data - trend

SA <- auto_SA(y, dates, detrend = FALSE)

y <- ts(data, start = c(1970,1), frequency = 12)
aa <- auto.arima(y)

summary(aa) #0.423
SA$MSE #0.302



SA <- auto_SA(y, dates)

ts.plot(cbind(y, SA$data_sa), col = c("blue", "red"))


#
# Seas <- seas(y, x11 = "", transform.function = "none", x11.appendfcst = "no")
# plot(dates, (SA$data_sa), type = "l", lwd = 1, col = "red",
#      xlab = "Dates", ylab = "IP")
# lines(dates, (final(Seas)), lwd = 1, col = "blue")
#
# ts.plot(final(Seas))







#
# bob <- Arima(y, order = c(1,0,1), seasonal = c(1,0,1))
#
#
#
#
# ts.plot(data_monthly)
#
# summary(sa)
#
# Seas <- seas(y, x11 = "", transform.function = "none", x11.appendfcst = "no")
#
# sa_factor <-  series(Seas, "x11.adjustfac")
#
# ts.plot(cbind(final(Seas),est$Y_sa), col = c("blue", "red") )
#
# ts.plot(cbind(sue,est$Y_sa), col = c("blue", "red") )
#
# ts.plot(cbind(final(Seas),y_sa), col = c("blue", "red") )
#
#
#
#
# summary(Seas)
#
#
# ts.plot(cbind(sa_factor, SA$adj_factor), col = c("blue", "red"))
#
# bob <- log(y) - sa_factor
# sue <- log(y) - SA$adj_factor
#
# ts.plot(cbind(tail(bob, 20),tail(sue, 20)), col = c("blue", "red"))
#
#
# summary(Seas)
#
#
#
# dt <- initial_claims
#
# monthly <- agg_to_monthly(dt)
#
# data_monthly <- monthly$ICNSA
# dates_monthly <- as.Date(monthly$DATE)
#
# SA <- auto_SA(data = data_monthly, dates = dates_monthly)
#
# # (1) Do factors add up?
# zero <- data_monthly - SA$data_sa - SA$adj_factor
# max(zero) #yes
#
# #(2) test against other packages (other stuff works for monthly data)
# ##### test against the forecast package ######################
#
# library(forecast)
# data_ts <- ts(data_monthly, start = c(1967,1), frequency = 12)
#
# decomposed <- decompose(data_ts)
# sa <- forecast::seasadj(decomposed)
#
# #Compare results
#
# ts.plot(cbind(sa, SA$data_sa), col = c("red", "blue")) #no contest
#
# ##### test against the seasonal package ######################
#
# library(seasonal)
# sa <- seas(data_ts, x11 = "", transform.function = "none") #, x11.appendfcst = "yes"
# sa_factor <- series(sa, "x11.adjustfac") #test the factors
# zero <- data_ts - final(sa) - sa_factor
# max(zero) #check
# ts.plot(cbind(final(sa), SA$data_sa), col = c("red", "blue")) #very similar... can't say who's best by just looking
# #Which series offers better predictive power?
# first <- auto.arima(SA$data_sa[-(1:13)])
# second <- auto.arima(final(sa)[-(1:13)])
# print(first)
# print(second)
#
# #seasonal does better, but at least we're close!
#
#
#
