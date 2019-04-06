#' Automatically seasonally adjust high frequency or irregular frequency data
#'
#' @param data data in matrix format with time in rows
#' @param dates dates corresponding to data points if data is not ts_boxable()
#' @param forecast T/F return forecasted adjustment factors
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats loess median model.matrix na.exclude optim predict var
#' @importFrom utils head tail
#' @useDynLib seasonality
auto_SA <- function(data, dates = NULL, forecast = FALSE){

  if(is.null(dates)){
    stopifnot(requireNamespace("tsbox"))
    if(ts_boxable(data)){
      dates <- ts_df(data)$time
    }else{
      stop("If 'dates' is NULL data must be ts_boxable()")
    }
  }else{
      if(length(dates) != length(data)){
        stop("Length of date vector must correspond to columns in data")
      }
  }
  dates   <- as.Date(dates)
  data <- c(unclass(data))

  OUT <- auto_SA_core(data, dates, forecast)

  return(OUT)

}


arma_loss <- function(x, y, order, ar_lags, ma_lags){
  idx <- order[1]
  if(order[1]>0){
    p <- t(as.matrix(x[1:idx]))
  }else{
    p <- matrix(0,0,0)
  }
  idx <- idx+1
  if(order[2]>0){
    q <-  t(as.matrix(x[idx:(idx+order[2]-1)]))
  }else{
    q <- matrix(0,0,0)
  }
  idx <- idx+order[2]
  if(!is.null(ar_lags)){
    P <-  t(as.matrix(x[idx:(idx+NCOL(ar_lags)-1)]))
    ar_lags <- as.matrix(ar_lags)
    idx <- idx + NCOL(ar_lags)
  }else{
    P <- matrix(0,0,0)
    ar_lags <- matrix(0,0,0)
  }
  if(!is.null(ma_lags)){
    Q <-  t(as.matrix(x[idx:(idx+NCOL(ma_lags)-1)]))
    ma_lags <- as.matrix(ma_lags)
  }else{
    Q <- matrix(0,0,0)
    ma_lags <- matrix(0,0,0)
  }
  est <- SARMA(Y = y, p = p, q = q, P = P, Q = Q, P_lag = ar_lags, Q_lag = ma_lags)
  return(est$MSE)
}


arma_est <- function(x, y, order, ar_lags, ma_lags){
  idx <- order[1]
  if(order[1]>0){
    p <- t(as.matrix(x[1:idx]))
  }else{
    p <- matrix(0,0,0)
  }
  idx <- idx+1
  if(order[2]>0){
    q <-  t(as.matrix(x[idx:(idx+order[2]-1)]))
  }else{
    q <- matrix(0,0,0)
  }
  idx <- idx+order[2]
  if(!is.null(ar_lags)){
    P <-  t(as.matrix(x[idx:(idx+NCOL(ar_lags)-1)]))
    ar_lags <- as.matrix(ar_lags)
    idx <- idx + NCOL(ar_lags)
  }else{
    P <- matrix(0,0,0)
    ar_lags <- matrix(0,0,0)
  }
  if(!is.null(ma_lags)){
    Q <-  t(as.matrix(x[idx:(idx+NCOL(ma_lags)-1)]))
    ma_lags <- as.matrix(ma_lags)
  }else{
    Q <- matrix(0,0,0)
    ma_lags <- matrix(0,0,0)
  }
  est <- SARMA(Y = y, p = p, q = q, P = P, Q = Q, P_lag = ar_lags, Q_lag = ma_lags)
  return(est)
}

make_lag_index <- function(dates, effects){
  if("year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = 1))
    lags <- c(which_date_closest_ordered(FromVec = year_ago, IndVec = dates)) - 1 #C++ indexing
  }else if("two_year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = 2))
    lags <- c(which_date_closest_ordered(FromVec = year_ago, IndVec = dates)) - 1 #C++ indexing
  }else if("three_year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = 3))
    lags <- c(which_date_closest_ordered(FromVec = year_ago, IndVec = dates)) - 1 #C++ indexing
  }else{
    lags <- NULL
  }
  return(lags)
}


seasonal_arma <- function(data, dates, fc_dates = NULL, order = c(1,1), seasonal_ar = NULL, seasonal_ma = NULL){

  ddates <- median(diff(dates))
  if(ddates>66){
    stop("Data must have a frequency between daily and quarterly")
  }
  if(tail(dates,1) - dates[1] < 1095){
    stop("At least three years of data are required for seasonal adjustment")
  }

  p <- matrix(0,1,order[1])
  q <- matrix(0,1,order[2])

  all_dates <- c(dates, fc_dates)

  if(!is.null(seasonal_ar)){
    if(is.numeric(seasonal_ar)){
      ar_lags <- seq(length(all_dates))%x%matrix(1,1,length(seasonal_ar)) - matrix(1,length(all_dates),1)%x%matrix(seasonal_ar, 1, length(seasonal_ar))
      ar_lags <- ar_lags-1 #c++ indexing
      ar_lags(ar_lags<0) <- 0
    }else{
      ar_lags <- do.call("cbind",lapply(seasonal_ar, FUN = make_lag_index, dates = all_dates))
    }
    P <- matrix(0,1,NCOL(ar_lags))
  }else{
    ar_lags <- NULL
    P <- NULL
  }
  if(!is.null(seasonal_ma)){
    if(is.numeric(seasonal_ma)){
      ma_lags <- seq(length(all_dates))%x%matrix(1,1,length(seasonal_ma)) - matrix(1,length(all_dates),1)%x%matrix(seasonal_ma, 1, length(seasonal_ma))
      ma_lags <- ar_lags-1 #c++ indexing
      ma_lags(ar_lags<0) <- 0
    }else{
      ma_lags <- do.call("cbind",lapply(seasonal_ma, FUN = make_lag_index, dates = all_dates))
    }
    Q <- matrix(0,1,NCOL(ma_lags))
  }else{
    ma_lags <-NULL
    Q <- NULL
  }

  x <- c(p,q,P,Q)
  if(length(x)==1){
    params <- optim(par=x, fn = arma_loss, y=data, order=order, ar_lags=ar_lags, ma_lags=ma_lags, method = "Brent", lower = -1, upper = 1)
  }else{
    params <- optim(par=x, fn = arma_loss, y=data, order=order, ar_lags=ar_lags, ma_lags=ma_lags)
  }
  est <- arma_est(x=params$par, y=data, order=order, ar_lags=ar_lags, ma_lags=ma_lags)

  #Format output
  est$Y_sa <- c(est$Y_sa)
  names(est$Y_sa) <- dates

  est$seas <- c(est$seas)
  names(est$seas) <- all_dates

  est$E <- c(est$E)
  names(est$E) <- dates

  est$YP <- c(est$YP)
  names(est$YP) <- dates

  return(est)
}


select_SARMA <- function(y,dates){
  ar_effect <- NULL #need to add holidays and such
  ma_effect <- NULL
  ar_it <- 1
  ma_it <- 1
  order <- c(1,0)
  order_new <- c(1,0)
  years <- c("year", "two_year", "three_year")
  est <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = ma_effect)
  ucv <- var(y, na.rm = T)
  pen <- .2
  SIC0 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect))/sqrt(length(y)) #Seth's bogus information criteria
  for(it in 1:4){
    ar_effect_new <- c(ar_effect, years[ar_it])
    est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect_new, seasonal_ma = ma_effect)
    SIC1 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect))/sqrt(length(y))  #Seth's bogus information criteria

    ma_effect_new <- c(ma_effect, years[ma_it])
    est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = ma_effect_new)
    SIC2 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect_new))/sqrt(length(y))  #Seth's bogus information criteria

    order_ar <- order + c(1,0)
    est <- seasonal_arma(y, dates, order = order_ar, seasonal_ar = ar_effect, seasonal_ma = ma_effect)
    SIC3 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect))/sqrt(length(y))  #Seth's bogus information criteria

    if(length(ma_effect)==0){ #q>0 requires a seasonal MA component
      order_ma <- order + c(0,1)
      ma_effect_new <- c(ma_effect, years[1])
      est <- seasonal_arma(y, dates, order = order_ma, seasonal_ar = ar_effect, seasonal_ma = ma_effect_new)
      SIC4 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect))/sqrt(length(y))  #Seth's bogus information criteria
      if( all(SIC0 > c(SIC1, SIC2, SIC3, SIC4)) ){
        break
      }else if( all(SIC1 > c(SIC0, SIC2, SIC3, SIC4)) ){
        ar_effect <- ar_effect_new
        ar_it <- ar_it + 1
        SIC0 <- SIC1
      }else if( all(SIC2 > c(SIC1, SIC0, SIC3, SIC4)) ){
        ma_effect <- ma_effect_new
        ma_it <- ma_it + 1
        SIC0 <- SIC2
      }else if( all(SIC3 > c(SIC1, SIC2, SIC0, SIC4)) ){
        order <- order_ar
        SIC0 <- SIC3
      }
      else if( all(SIC4 > c(SIC0, SIC1, SIC2, SIC3)) ){
        order <- order_ma
        ma_effect <- ma_effect_new
        ma_it <- 2
        SIC0 <- SIC4
      }
    }else{ #if length(ma_effect>0)
      order_ma <- order + c(0,1)
      est <- seasonal_arma(y, dates, order = order_ma, seasonal_ar = ar_effect, seasonal_ma = ma_effect)
      SIC4 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect))/sqrt(length(y))  #Seth's bogus information criteria
      if( all(SIC0 > c(SIC1, SIC2, SIC3, SIC4)) ){
        break
      }else if( all(SIC1 > c(SIC0, SIC2, SIC3, SIC4)) ){
        ar_effect <- ar_effect_new
        ar_it <- ar_it + 1
        SIC0 <- SIC1
      }else if( all(SIC2 > c(SIC1, SIC0, SIC3, SIC4)) ){
        ma_effect <- ma_effect_new
        ma_it <- ma_it + 1
        SIC0 <- SIC2
      }else if( all(SIC3 > c(SIC1, SIC2, SIC0, SIC4)) ){
        order <- order_ar
        SIC0 <- SIC3
      }
      else if( all(SIC4 > c(SIC0, SIC1, SIC2, SIC3)) ){
        order <- order_ma
        SIC0 <- SIC4
      }
    }#if length(ma_effect>0)
  } #for(it in 1:4)

  out <- list(
    order = order,
    ar_effect = ar_effect,
    ma_effect = ma_effect
  )
  return(out)
}

auto_SA_core <- function(data, dates, forecast = FALSE){

  take_logs <- all(data>0)
  if(take_logs){
    data <- log(data)
  }

  dates <- as.Date(dates)
  ddates <- median(diff.Date(dates))
  if(ddates >= 28 && ddates <= 31){ #Monthly data
    trend <- loess(data ~ seq(length(data)), span = .12, na.action = na.exclude) #refine parameter selection
    trend <- predict(trend)
    #ts.plot(cbind(trend, data), col = c("red", "blue"))
    y <- data - trend
    dates <- first_of_month(dates)
    if(is.numeric(forecast)){ #forecast must be numeric or logical
      fc_dates = forecast
    }else if(is.logical(forecast)){
      if(forecast){
        fc_dates <- seq.Date(from = tail(dates,1), by = "month", length.out = 13)[-1]
      }else{
        fc_dates <- NULL
      }
    }else{
      fc_dates = NULL
    }
    #Step (1): select SARMA model
    params <- select_SARMA(y,dates)
    est <- seasonal_arma(data = y, dates = dates, fc_dates = fc_dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect)

    #Step (2): Select holiday effects (to be implemented!)

    #Step (3): Remove weeday/trading day effects and other exogenous junk (just trading days so far)
    trading_days <- do.call("c",lapply(dates, FUN = weekdays_in_month))
    td <- UVreg(as.matrix(trading_days), y = est$E, rm_outlier = 1)
    if(abs(td$B/td$sd)>2){
      trading_days_long <- do.call("c",lapply(c(dates, fc_dates), FUN = weekdays_in_month))
      est$seas <- est$seas + td$B*trading_days_long
      est$Y_sa <- est$Y_sa - td$B*trading_days
      est$trading_days <- TRUE
    }
  } #If data is monthly

  if(ddates >= 6 && ddates <= 8){ #Weekly data
    trend <- loess(data ~ seq(length(data)), span = .15, na.action = na.exclude) #refine parameter selection
    trend <- predict(trend) #all we're really interested in
    #ts.plot(cbind(data,trend), col = c("red", "blue"))
    y <- data - trend

    if(is.logical(forecast)){
      if(forecast){
        if(all(unique(ddates)==7)){ #true weekly
          fc_dates <- seq.Date(from = tail(dates,1), by = 7, length.out = 25)[-1]
        }else if(all(sapply(dates, FUN = day)%in%c(7,14,21,28:31))){ #pseudo weekly
          fc_dates <- pseudo_weekly_sequence(start = tail(dates,1), length = 25)
        }else{
          fc_dates <- NULL
        }
      }else{
        fc_dates <- NULL
      }
    }else{
      fc_dates = NULL
    }

    #Step (1): select SARMA model
    params <- select_SARMA(y,dates)
    est <- seasonal_arma(data = y, dates = dates, fc_dates = fc_dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect)

    #Step (2): Select holiday effects (to be implemented!)

    #Step (3): Remove weeday/trading day effects and other exogenous junk (just trading days so far)

  } #if data is weekly

  if(ddates == 1){ #daily data
    trend <- loess(data ~ seq(length(data)), span = .3, na.action = na.exclude) #refine parameter selection
    trend <- predict(trend) #all we're really interested in
    #ts.plot(cbind(data,trend), col = c("red", "blue"))
    y <- data - trend

    if(is.numeric(forecast)){ #forecast must be numeric or logical
      fc_dates = forecast
    }else if(is.logical(forecast)){
      if(forecast){
        fc_dates <- seq.Date(from = tail(dates,1), by = 1, length.out = 61)[-1]
        }
    }else{
      fc_dates = NULL
    }

    #Step (1): select SARMA model
    params <- select_SARMA(y,dates)
    est <- seasonal_arma(data = y, dates = dates, fc_dates = fc_dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect)

    #Step (2): Select holiday effects (to be implemented!)

    #Step (3): Remove weeday/trading day effects and other exogenous junk (just trading days so far)

  } #if data is daily

  data_sa = data - est$seas[1:length(y)]
  if(take_logs){
    data_sa = exp(data_sa)
  }

  OUT <- list(
    data_sa = data_sa,
    adj_factor = est$seas,
    parameters = list(
      p = est$p,
      q = est$q,
      P = est$P,
      Q = est$Q
    ),
    log_transform = take_logs
  )
  return(OUT)

}


