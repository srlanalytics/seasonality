dummy_matrix <- function(x) {
  # factors with same level order as x
  xf <- factor(x, levels = unique(x))
  z <- model.matrix(~ xf + 0)
  attr(z, "assign") <- NULL
  attr(z, "contrasts") <- NULL
  colnames(z) <- levels(xf)
  # z[,-1]
  z
}

genhol_dates <- function(x, dates, start = 0, end = 0) {
  stopifnot(inherits(x, "Date"))
  stopifnot(inherits(dates, "Date"))
  freq <- which_freq(dates)
  timestamp <- timestamp_fun(freq)
  sq <- daily_seq(dates)


  idx.holiday <- which(sq %in% x)
  st <- idx.holiday + start
  en <- idx.holiday + end
  idx.holiday.extended <- unlist(Map(function(st, en) st:en, st, en))

  is.holiday <- logical(length(sq))
  is.holiday[idx.holiday.extended] <- TRUE

  z <- as.matrix(tapply(is.holiday, timestamp(sq), sum))
  colnames(z) <- "holiday"
  z
}

weekdays_in_month <- function(date){
  sq <- seq.Date(from = first_of_month(date), to = end_of_month(date), by = "day")
  count <- sum(as.integer(!(weekdays(sq) %in% c("Saturday", "Sunday"))))
  if(as.numeric(format.Date(date,"%m"))==12){
    count = count - sum(as.integer(!(weekdays(sq[24:25]) %in% c("Saturday", "Sunday")))) #subtract Christmas eve and Christmas if counted
  }
  if(as.numeric(format.Date(date,"%m"))==1){
    count = count - as.integer(!(weekdays(sq[1]) %in% c("Saturday", "Sunday"))) #subtract New Years if counted
  }
  return(count)
}

closest_weekday <- function(dte, weekday = "Sunday"){
  sq   <- c(seq.Date(from = dte-6, by = 1, length.out = 6), seq.Date(from = dte, by = 1, length.out = 7))
  pos  <- (which(weekdays(sq) == weekday) - 7) #position of Fridays relative to the current date
  date_out <- sq[ pos[which.min(abs(pos))]+7 ]
  return(date_out)
}

previous_weekday <- function(dte, weekday = "Sunday"){
  sq   <- seq.Date(from = dte-6, by = 1, length.out = 7)
  date_out <- sq[which(weekdays(sq) == weekday)]
  return(date_out)
}

# auto detect frequency
get_freq_dates <- function(y, dates) {
  out <- median(diff(dates[which(!is.na(y))]))
}

# auto detect frequency
get_freq <- function(dates) {
  out <- median(diff(which(!is.na(dates))))
}
