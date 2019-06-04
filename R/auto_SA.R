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
auto_SA <- function(data, dates = NULL, forecast = FALSE, take_logs = "auto"){

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


#which(dates == previous_year[this_year == dates[993]])

seasonal_rev_arma <- function(data, dates, order = c(1,0), seasonal_ar = NULL, seasonal_ma = NULL, initial_vals = NULL){

  ddates <- median(diff(dates))
  if(ddates>66){
    stop("Data must have a frequency between daily and quarterly")
  }
  if(tail(dates,1) - dates[1] < 1095){
    stop("At least three years of data are required for seasonal adjustment")
  }

  if(is.null(initial_vals)){
    p <- matrix(0,1,order[1])
    q <- matrix(0,1,order[2])
  }else{
    p <- initial_vals$p
    q <- initial_vals$q
  }

  if(!is.null(seasonal_ar)){
    ar_lags <- make_lead_index(dates, effects = seasonal_ar)
    if(is.null(initial_vals)){
      P <- matrix(0,1,NCOL(ar_lags))
    }else{
      P <- initial_vals$P
    }
  }else{
    ar_lags <- NULL
    P <- NULL
  }
  if(!is.null(seasonal_ma)){
    ma_lags <- make_lag_index(dates, effects = seasonal_ma)
    if(is.null(initial_vals)){
      Q <- matrix(0,1,NCOL(ma_lags))
    }else{
      Q <- initial_vals$Q
    }
  }else{
    ma_lags <-NULL
    Q <- NULL
  }

  x <- c(p,q,P,Q)
  if(length(x)==1){
    params <- optim(par=x, fn = arma_loss, y=rev(data), order=order, ar_lags=ar_lags, ma_lags=ma_lags, method = "Brent", lower = -1, upper = 1)
  }else{
    params <- optim(par=x, fn = arma_loss, y=rev(data), order=order, ar_lags=ar_lags, ma_lags=ma_lags)
  }
  est <- arma_est(x=params$par, y=rev(data), order=order, ar_lags=ar_lags, ma_lags=ma_lags)

  #Format output
  est$Y_sa <- rev(c(est$Y_sa))
  names(est$Y_sa) <- dates

  if(is.null(ar_lags) && is.null(ma_lags)){
    est$seas <- rep(0,length(dates))
    names(est$seas) <- dates
  }else{
    est$seas <- rev(c(est$seas))
    names(est$seas) <- dates
  }

  est$E <- rev(c(est$E))
  names(est$E) <- dates

  est$YP <- rev(c(est$YP))
  names(est$YP) <- dates

  return(est)
}

seasonal_arma <- function(data, dates, fc_dates = NULL, order = c(1,0), seasonal_ar = NULL, seasonal_ma = NULL, initial_vals = NULL){

  ddates <- median(diff(dates))
  if(ddates>66){
    stop("Data must have a frequency between daily and quarterly")
  }
  if(tail(dates,1) - dates[1] < 1095){
    stop("At least three years of data are required for seasonal adjustment")
  }

  if(is.null(initial_vals)){
    p <- matrix(0,1,order[1])
    q <- matrix(0,1,order[2])
  }else{
    p <- initial_vals$p
    q <- initial_vals$q
  }

  all_dates <- c(dates, fc_dates)

  if(!is.null(seasonal_ar)){
    if(is.numeric(seasonal_ar)){
      ar_lags <- seq(length(all_dates))%x%matrix(1,1,length(seasonal_ar)) - matrix(1,length(all_dates),1)%x%matrix(seasonal_ar, 1, length(seasonal_ar))
      ar_lags <- ar_lags-1 #c++ indexing
      ar_lags(ar_lags<0) <- 0
    }else{
      ar_lags <- make_lag_index(all_dates, effects = seasonal_ar)
    }
    if(is.null(initial_vals)){
      P <- matrix(0,1,NCOL(ar_lags))
    }else{
      P <- initial_vals$P
    }
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
      ma_lags <- make_lag_index(all_dates, effects = seasonal_ma)
    }
    if(is.null(initial_vals)){
      Q <- matrix(0,1,NCOL(ma_lags))
    }else{
      Q <- initial_vals$Q
    }
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

  if(is.null(ar_lags) && is.null(ma_lags)){
    est$seas <- rep(0,length(all_dates))
    names(est$seas) <- all_dates
  }else{
    est$seas <- c(est$seas)
    names(est$seas) <- all_dates
  }

  est$E <- c(est$E)
  names(est$E) <- dates

  est$YP <- c(est$YP)
  names(est$YP) <- dates

  return(est)
}

auto_SA_core <- function(data, dates, forecast = FALSE, take_logs = "auto"){

  if(take_logs == "auto"){
    take_logs <- all(data[is.finite(data)]>0)
  }
  if(take_logs){
    data <- log(data)
  }

  dates <- as.Date(dates)
  ddates <- median(diff.Date(dates))
  if(ddates >= 28 && ddates <= 31){ #Monthly data
    trend <- loess(data ~ seq(length(data)), span = .1, na.action = na.exclude) #refine parameter selection
    trend <- predict(trend)
    #ts.plot(cbind(trend, data), col = c("red", "blue"))
    y <- data - trend
    dates <- end_of_month(dates)
    if(is.numeric(forecast)){ #forecast must be numeric or logical
      fc_dates <- end_of_month(seq.Date(from = first_of_month(tail(dates,1)), by = "month", length.out = forecast)[-1])
    }else if(is.logical(forecast)){
      if(forecast){
        fc_dates <- end_of_month(seq.Date(from = first_of_month(tail(dates,1)), by = "month", length.out = 13)[-1])
      }else{
        fc_dates <- NULL
      }
    }else{
      fc_dates = NULL
    }
    #Step (1): select SARMA model
    params <- select_SARMA(y,dates,fc_dates)
    #est <- seasonal_arma(data = y, dates = dates, fc_dates = fc_dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect)

    #Step (2): Smooth --- need to do the properly at some point
    #Question: will forecasts made on smoothed data perform better since the tail will not be smoothed (i.e. is just filtered)?
    initial_vals <- list(p = params$p,
                         q = params$q,
                         P = params$P,
                         Q = params$Q)
    rev_est <- seasonal_rev_arma(y, dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect,
                                 initial_vals = initial_vals)
    idx <- seq(1,length(dates)) #if fc_dates, params$seas will be longer than y
    sa <- params$seas[idx]
    sa[sa==0] <- rev_est$seas[sa==0]
    params$seas[idx] <- sa
    y_sa <- y - sa
    cp <- ncol(params$p)
    z <- stack_obs(as.matrix(y_sa), cp)
    E <- c(rep(0, cp), y_sa[-seq(1,cp)] - params$p%*%z[-nrow(z),])

    #Step (3): Remove weekday/trading day effects and other exogenous junk (just trading days so far)

    #trading_days

    trading_days <- do.call("c",lapply(dates, FUN = weekdays_in_month)) #need to update to inlcude weekdays in quarter for quarterly data
    trading_days <- scale(trading_days)
    bet <- mean(trading_days*E)
    Enew <- E - bet*trading_days
    if( (mean(E^2) - mean(Enew^2))/var(y) > .1/sqrt(length(y))){
      td <- do.call("c",lapply(as.Date(names(params$seas)), FUN = weekdays_in_month)) #this clunky approach accomadates fc_dates
      td <- (td - attr(trading_days, "scaled:center"))/attr(trading_days, "scaled:scale")
      params$seas <- params$seas + bet*td
      y_sa <- y_sa - bet*trading_days
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
    params <- select_SARMA(y,dates, fc_dates)
    #est <- seasonal_arma(data = y, dates = dates, fc_dates = fc_dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect)

    #Step (2): Smooth --- fudged for now...
    initial_vals <- list(p = params$p,
                         q = params$q,
                         P = params$P,
                         Q = params$Q)
    rev_est <- seasonal_rev_arma(y, dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect,
                                 initial_vals = initial_vals)
    idx <- seq(1,length(dates)) #if fc_dates, params$seas will be longer than y
    sa <- params$seas[idx]
    sa[sa==0] <- rev_est$seas[sa==0]
    params$seas[idx] <- sa
    # y_sa <- y - sa
    # cp <- ncol(params$p)
    # z <- stack_obs(as.matrix(y_sa), cp)
    # E <- c(rep(0, cp), y_sa[-seq(1,cp)] - params$p%*%z[-nrow(z),])

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
    params <- select_SARMA(y,dates, fc_dates)
    #est <- seasonal_arma(data = y, dates = dates, fc_dates = fc_dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect)

    #Step (2): Smooth... fudged
    initial_vals <- list(p = params$p,
                         q = params$q,
                         P = params$P,
                         Q = params$Q)
    rev_est <- seasonal_rev_arma(y, dates, order = params$order, seasonal_ar = params$ar_effect, seasonal_ma = params$ma_effect,
                                 initial_vals = initial_vals)
    idx <- seq(1,length(dates)) #if fc_dates, params$seas will be longer than y
    sa <- params$seas[idx]
    sa[sa==0] <- rev_est$seas[sa==0]
    params$seas[idx] <- sa
    # y_sa <- y - sa
    # cp <- ncol(params$p)
    # z <- stack_obs(as.matrix(y_sa), cp)
    # E <- c(rep(0, cp), y_sa[-seq(1,cp)] - params$p%*%z[-nrow(z),])

    #Step (3): Remove weeday/trading day effects and other exogenous junk (just trading days so far)

  } #if data is daily

  data_sa <- data - sa

  if(take_logs){
    data_sa = exp(data_sa)
  }

  OUT <- list(
    data_sa = data_sa,
    adj_factor = params$seas,
    parameters = list(
      p = params$p,
      q = params$q,
      P = params$P,
      Q = params$Q
    ),
    log_transform = take_logs
  )
  return(OUT)

}


