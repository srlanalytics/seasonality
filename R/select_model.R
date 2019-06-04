
holiday_range <- function(hdate, dates, before = 30, after = 30){
  all_dates <- seq.Date(from = hdate - before, to = hdate + after, by = "day")
  names(all_dates) <- seq(-before,after)
  out_dates <- all_dates[all_dates%in%dates]
  return(out_dates)
}

min_diff_name <- function(this_date_name, last_year_names) which.min(abs(as.numeric(last_year_names) - as.numeric(this_date_name) ))

previous_year_holiday <- function(idx, Hdates){
  this_year <- Hdates[[idx]]
  last_year <- Hdates[[idx-1]]
  out <- last_year[sapply(names(this_year), FUN = min_diff_name, last_year_names = names(last_year))]
  return(out)
}

which_is_equal <- function(x,y) which(y==x) #x - single value, y - vector to compare to

make_lead_index <- function(dates, effects){
  lags <- NULL
  if("year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = -1))
    tmp  <- rev(length(dates) - which_date_closest_ordered(FromVec = year_ago, IndVec = dates) + 1)
    tmp[seq(1,max(which(tmp==1))-1)] <- 0
    lags <- cbind(lags, tmp) #R indexing
  }
  if("two_year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = -2))
    tmp  <- rev(length(dates) - which_date_closest_ordered(FromVec = year_ago, IndVec = dates) + 1)
    tmp[seq(1,max(which(tmp==1))-1)] <- 0
    lags <- cbind(lags, tmp) #R indexing
  }
  if("three_year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = -3))
    tmp  <- rev(length(dates) - which_date_closest_ordered(FromVec = year_ago, IndVec = dates) + 1)
    tmp[seq(1,max(which(tmp==1))-1)] <- 0
    lags <- cbind(lags, tmp) #R indexing
  }
  dates <- rev(dates) #above dates must be ordered for which_date_closest_ordered function. Here that's not neccessary.
  if("cny"%in%effects){
    tmp <- rep(0, length(dates))
    hdates <- rev(cny[cny>= min(dates) & cny<=max(dates)])
    Hdates <- lapply(hdates, FUN = holiday_range, dates = dates)
    previous_year <- do.call("c",  lapply(seq(2,length(Hdates)), FUN = previous_year_holiday, Hdates = Hdates))
    this_year <-  do.call("c", Hdates)[-seq(1,length(Hdates[[1]]))]
    #if(is.unsorted(dates, na.rm = TRUE) || is.unsorted(this_year, na.rm = TRUE) || is.unsorted(previous_year, na.rm = TRUE)) stop("Dates must be sorted")
    tmp[sapply(this_year, FUN = which_is_equal, y = dates)] <- sapply(previous_year, FUN = which_is_equal, y = dates)
    lags <- cbind(lags, tmp) #R indexing
  }
  if("easter"%in%effects){
    tmp <- rep(0, length(dates))
    hdates <- rev(easter[easter>= min(dates) & easter<=max(dates)])
    if(abs(median(diff(dates)))==31){
      Hdates <- lapply(hdates, FUN = holiday_range, dates = dates, before = 30, after = 30)
    }else{
      Hdates <- lapply(hdates, FUN = holiday_range, dates = dates, before = 7, after = 7)
    }
    previous_year <- do.call("c",  lapply(seq(2,length(Hdates)), FUN = previous_year_holiday, Hdates = Hdates))
    this_year <-  do.call("c", Hdates)[-seq(1,length(Hdates[[1]]))]
    #if(is.unsorted(dates, na.rm = TRUE) || is.unsorted(this_year, na.rm = TRUE) || is.unsorted(previous_year, na.rm = TRUE)) stop("Dates must be sorted")
    tmp[sapply(this_year, FUN = which_is_equal, y = dates)] <- sapply(previous_year, FUN = which_is_equal, y = dates)
    lags <- cbind(lags, tmp) #R indexing
  }
  if("diwali"%in%effects){
    tmp <- rep(0, length(dates))
    hdates <- rev(diwali[diwali>= min(dates) & diwali<=max(dates)])
    Hdates <- lapply(hdates, FUN = holiday_range, dates = dates)
    previous_year <- do.call("c",  lapply(seq(2,length(Hdates)), FUN = previous_year_holiday, Hdates = Hdates))
    this_year <-  do.call("c", Hdates)[-seq(1,length(Hdates[[1]]))]
    #if(is.unsorted(dates, na.rm = TRUE) || is.unsorted(this_year, na.rm = TRUE) || is.unsorted(previous_year, na.rm = TRUE)) stop("Dates must be sorted")
    tmp[sapply(this_year, FUN = which_is_equal, y = dates)] <- sapply(previous_year, FUN = which_is_equal, y = dates)
    lags <- cbind(lags, tmp) #R indexing
  }
  if("black_friday"%in%effects){
    tmp <- rep(0, length(dates))
    hdates <- rev(black_friday[black_friday>= min(dates) & black_friday<=max(dates)])
    Hdates <- lapply(hdates, FUN = holiday_range, dates = dates, before = 6, after = 6)
    previous_year <- do.call("c",  lapply(seq(2,length(Hdates)), FUN = previous_year_holiday, Hdates = Hdates))
    this_year <-  do.call("c", Hdates)[-seq(1,length(Hdates[[1]]))]
    #if(is.unsorted(dates, na.rm = TRUE) || is.unsorted(this_year, na.rm = TRUE) || is.unsorted(previous_year, na.rm = TRUE)) stop("Dates must be sorted")
    tmp[sapply(this_year, FUN = which_is_equal, y = dates)] <- sapply(previous_year, FUN = which_is_equal, y = dates)
    lags <- cbind(lags, tmp) #R indexing
  }

  return(lags) #keep R indexing
}

make_lag_index <- function(dates, effects){
  lags <- NULL
  if("year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = 1))
    tmp  <- which_date_closest_ordered(FromVec = year_ago, IndVec = dates)
    tmp[seq(1,max(which(tmp==1))-1)] <- 0
    lags <- cbind(lags, tmp) #R indexing
  }
  if("two_year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = 2))
    tmp <- which_date_closest_ordered(FromVec = year_ago, IndVec = dates)
    tmp[seq(1,max(which(tmp==1))-1)] <- 0
    lags <- cbind(lags, tmp) #R indexing
  }
  if("three_year"%in%effects){
    #this is for an irregular patern of dates
    year_ago <- do.call("c", lapply(dates, years_ago, shift = 3))
    tmp <- which_date_closest_ordered(FromVec = year_ago, IndVec = dates)
    tmp[seq(1,max(which(tmp==1))-1)] <- 0
    lags <- cbind(lags, tmp) #R indexing
  }
  if("cny"%in%effects){
    tmp <- rep(0, length(dates))
    hdates <- cny[cny>= min(dates) & cny<=max(dates)]
    Hdates <- lapply(hdates, FUN = holiday_range, dates = dates)
    previous_year <- do.call("c",  lapply(seq(2,length(Hdates)), FUN = previous_year_holiday, Hdates = Hdates))
    this_year <-  do.call("c", Hdates)[-seq(1,length(Hdates[[1]]))]
    #if(is.unsorted(dates, na.rm = TRUE) || is.unsorted(this_year, na.rm = TRUE) || is.unsorted(previous_year, na.rm = TRUE)) stop("Dates must be sorted")
    tmp[sapply(this_year, FUN = which_is_equal, y = dates)] <- sapply(previous_year, FUN = which_is_equal, y = dates)
    lags <- cbind(lags, tmp) #R indexing
  }
  if("easter"%in%effects){
    tmp <- rep(0, length(dates))
    hdates <- easter[easter>= min(dates) & easter<=max(dates)]
    if(median(diff(dates))==31){
      Hdates <- lapply(hdates, FUN = holiday_range, dates = dates, before = 30, after = 30)
    }else{
      Hdates <- lapply(hdates, FUN = holiday_range, dates = dates, before = 7, after = 7)
    }
    previous_year <- do.call("c",  lapply(seq(2,length(Hdates)), FUN = previous_year_holiday, Hdates = Hdates))
    this_year <-  do.call("c", Hdates)[-seq(1,length(Hdates[[1]]))]
    #if(is.unsorted(dates, na.rm = TRUE) || is.unsorted(this_year, na.rm = TRUE) || is.unsorted(previous_year, na.rm = TRUE)) stop("Dates must be sorted")
    tmp[sapply(this_year, FUN = which_is_equal, y = dates)] <- sapply(previous_year, FUN = which_is_equal, y = dates)
    lags <- cbind(lags, tmp) #R indexing
  }
  if("diwali"%in%effects){
    tmp <- rep(0, length(dates))
    hdates <- diwali[diwali>= min(dates) & diwali<=max(dates)]
    Hdates <- lapply(hdates, FUN = holiday_range, dates = dates)
    previous_year <- do.call("c",  lapply(seq(2,length(Hdates)), FUN = previous_year_holiday, Hdates = Hdates))
    this_year <-  do.call("c", Hdates)[-seq(1,length(Hdates[[1]]))]
    #if(is.unsorted(dates, na.rm = TRUE) || is.unsorted(this_year, na.rm = TRUE) || is.unsorted(previous_year, na.rm = TRUE)) stop("Dates must be sorted")
    tmp[sapply(this_year, FUN = which_is_equal, y = dates)] <- sapply(previous_year, FUN = which_is_equal, y = dates)
    lags <- cbind(lags, tmp) #R indexing
  }
  if("black_friday"%in%effects){
    tmp <- rep(0, length(dates))
    hdates <- black_friday[black_friday>= min(dates) & black_friday<=max(dates)]
    Hdates <- lapply(hdates, FUN = holiday_range, dates = dates, before = 6, after = 6)
    previous_year <- do.call("c",  lapply(seq(2,length(Hdates)), FUN = previous_year_holiday, Hdates = Hdates))
    this_year <-  do.call("c", Hdates)[-seq(1,length(Hdates[[1]]))]
    #if(is.unsorted(dates, na.rm = TRUE) || is.unsorted(this_year, na.rm = TRUE) || is.unsorted(previous_year, na.rm = TRUE)) stop("Dates must be sorted")
    tmp[sapply(this_year, FUN = which_is_equal, y = dates)] <- sapply(previous_year, FUN = which_is_equal, y = dates)
    lags <- cbind(lags, tmp) #R indexing
  }

  return(lags) #keep R indexing
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
  pen <- .15
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
  } #for(it in 1:4)

  #test holiday effects --- need to update to set initial values to results from previous iteration
  #cny

  est <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = ma_effect)
  #est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = ma_effect)
  # ar_lags <- make_lag_index(dates, effects = "year")
  # ma_lags <- make_lag_index(dates, effects = c("year", "cny"))
  # est <- SARMA(Y = y, p = est$p, q = est$q, P = est$P, Q = est$Q, P_lag = ar_lags, Q_lag = ma_lags)

  initial_vals <- list(p = est$p,
                       q = est$q,
                       P = est$P,
                       Q = est$Q)
  new_model <- FALSE
  ar_effect_new <- ar_effect
  ma_effect_new <- ma_effect
  initial_tmp <- initial_vals
  initial_tmp$P <- c(initial_tmp$P, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = c(ar_effect, "cny"), seasonal_ma = ma_effect,
                        initial_vals = initial_tmp)
  SIC1 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 1)/sqrt(length(y))
  initial_tmp <- initial_vals
  initial_tmp$Q <- c(initial_tmp$Q, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = c(ma_effect, "cny"),
                        initial_vals = initial_tmp)
  SIC2 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 1)/sqrt(length(y))
  initial_tmp$P <- c(initial_tmp$P, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = c(ma_effect, "cny"), seasonal_ma = c(ma_effect, "cny"),
                        initial_vals = initial_tmp)
  SIC3 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 2)/sqrt(length(y))

  if(SIC1>SIC0 && SIC1>SIC2 && SIC1>SIC3){
    ar_effect_new <- c(ar_effect_new, "cny")
    SIC0 <- SIC1
    new_model <- TRUE
  }else if(SIC2>SIC0 && SIC2>SIC1 && SIC2>SIC3){
    ma_effect_new <- c(ma_effect_new, "cny")
    SIC0 <- SIC2
    new_model <- TRUE
  }else if(SIC3>SIC0 && SIC3>SIC2 && SIC3>SIC1){
    ar_effect_new <- c(ar_effect_new, "cny")
    ma_effect_new <- c(ma_effect_new, "cny")
    SIC0 <- SIC3
    new_model <- TRUE
  }

  #diwali

  initial_tmp <- initial_vals
  initial_tmp$P <- c(initial_tmp$P, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = c(ar_effect, "diwali"), seasonal_ma = ma_effect, initial_vals = initial_tmp)
  SIC1 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 1)/sqrt(length(y))
  initial_tmp <- initial_vals
  initial_tmp$Q <- c(initial_tmp$Q, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = c(ma_effect, "diwali"), initial_vals = initial_tmp)
  SIC2 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 1)/sqrt(length(y))
  initial_tmp$P <- c(initial_tmp$P, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = c(ma_effect, "diwali"), seasonal_ma = c(ma_effect, "diwali"), initial_vals = initial_tmp)
  SIC3 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 2)/sqrt(length(y))

  if(SIC1>SIC0 && SIC1>SIC2 && SIC1>SIC3){
    ar_effect_new <- c(ar_effect_new, "diwali")
    SIC0 <- SIC1
    new_model <- TRUE
  }else if(SIC2>SIC0 && SIC2>SIC1 && SIC2>SIC3){
    ma_effect_new <- c(ma_effect_new, "diwali")
    SIC0 <- SIC2
    new_model <- TRUE
  }else if(SIC3>SIC0 && SIC3>SIC2 && SIC3>SIC1){
    ar_effect_new <- c(ar_effect_new, "diwali")
    ma_effect_new <- c(ma_effect_new, "diwali")
    SIC0 <- SIC3
    new_model <- TRUE
  }

  #easter

  #est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = ma_effect)

  initial_tmp <- initial_vals
  initial_tmp$P <- c(initial_tmp$P, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = c(ar_effect, "easter"), seasonal_ma = ma_effect, initial_vals = initial_tmp)
  SIC1 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 1)/sqrt(length(y))
  initial_tmp <- initial_vals
  initial_tmp$Q <- c(initial_tmp$Q, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = c(ma_effect, "easter"), initial_vals = initial_tmp)
  SIC2 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 1)/sqrt(length(y))
  initial_tmp$P <- c(initial_tmp$P, 0)
  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = c(ma_effect, "easter"), seasonal_ma = c(ma_effect, "easter"), initial_vals = initial_tmp)
  SIC3 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 2)/sqrt(length(y))

  if(SIC1>SIC0 && SIC1>SIC2 && SIC1>SIC3){
    ar_effect_new <- c(ar_effect_new, "easter")
    SIC0 <- SIC1
    new_model <- TRUE
  }else if(SIC2>SIC0 && SIC2>SIC1 && SIC2>SIC3){
    ma_effect_new <- c(ma_effect_new, "easter")
    SIC0 <- SIC2
    new_model <- TRUE
  }else if(SIC3>SIC0 && SIC3>SIC2 && SIC3>SIC1){
    ar_effect_new <- c(ar_effect_new, "easter")
    ma_effect_new <- c(ma_effect_new, "easter")
    SIC0 <- SIC3
    new_model <- TRUE
  }


  #black_friday

  if(median(diff(dates))<8){
    est  <- seasonal_arma(y, dates, order = order, seasonal_ar = c(ar_effect, "black_friday"), seasonal_ma = ma_effect)
    SIC1 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 1)/sqrt(length(y))
    est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = c(ma_effect, "black_friday"))
    SIC2 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 1)/sqrt(length(y))
    est  <- seasonal_arma(y, dates, order = order, seasonal_ar = c(ma_effect, "black_friday"), seasonal_ma = c(ma_effect, "black_friday"))
    SIC3 <- (ucv-est$MSE)/ucv - pen*(sum(order) + length(ar_effect) + length(ma_effect) + 2)/sqrt(length(y))

    if(SIC1>SIC0 && SIC1>SIC2 && SIC1>SIC3){
      ar_effect_new <- c(ar_effect_new, "black_friday")
      SIC0 <- SIC1
      new_model <- TRUE
    }else if(SIC2>SIC0 && SIC2>SIC1 && SIC2>SIC3){
      ma_effect_new <- c(ma_effect_new, "black_friday")
      SIC0 <- SIC2
      new_model <- TRUE
    }else if(SIC3>SIC0 && SIC3>SIC2 && SIC3>SIC1){
      ar_effect_new <- c(ar_effect_new, "black_friday")
      ma_effect_new <- c(ma_effect_new, "black_friday")
      SIC0 <- SIC3
      new_model <- TRUE
    }

  }

  if(new_model){
    ar_effect <- ar_effect_new
    ma_effect <- ma_effect_new
  }

  est  <- seasonal_arma(y, dates, order = order, seasonal_ar = ar_effect, seasonal_ma = ma_effect)



  out <- list(
    order = order,
    ar_effect = ar_effect,
    ma_effect = ma_effect,
    p = est$p,
    q = est$q,
    P = est$P,
    Q = est$Q,
    seas = est$seas
    # Y_sa = est$Y_sa,
    # E = est$E,
    # YP = est$YP
  )
  return(out)
}