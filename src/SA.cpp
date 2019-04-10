// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

arma::sp_mat sp_rows(arma::sp_mat A,
                     arma::uvec r   ){
  uword n_rows   = A.n_rows;
  //  uword n_cols   = A.n_cols;
  uword n_r      = r.size();
  uvec  tmp      = regspace<uvec>(0,n_rows-1);
  tmp      = tmp.elem(r);
  umat  location = join_vert(trans(regspace<uvec>(0,n_r-1)),trans(tmp));
  sp_mat J(location,ones<vec>(n_r),n_r,n_rows);
  sp_mat C       = J*A;
  return(C);
}

arma::sp_mat sp_cols(arma::sp_mat A,
                     arma::uvec r   ){
  //  uword n_rows   = A.n_rows;
  uword n_cols   = A.n_cols;
  uword n_r      = r.size();
  uvec  tmp      = regspace<uvec>(0,n_cols-1);
  tmp            = tmp.elem(r);
  umat  location = join_vert(trans(tmp),trans(regspace<uvec>(0,n_r-1)));
  sp_mat J(location,ones<vec>(n_r),n_cols,n_r);
  sp_mat C       = A*J;
  return(C);
}

// [[Rcpp::export]]
arma::uvec get_ar_lags(arma::uword lag_length,
                       arma::uvec s_lag){
  uvec one(lag_length+1, fill::ones);
  uvec out = s_lag(0)*one - regspace<uvec>(0,lag_length);
  uvec ind;
  for(uword k=1; k<s_lag.size(); k++){
    out = join_vert(out, s_lag(k)*one - regspace<uvec>(0,lag_length));
  }
  return(out);
}

//Create the companion form of the transition matrix B
// [[Rcpp::export]]
arma::mat comp_form(arma::mat B){
  uword r = B.n_rows;
  uword c = B.n_cols;
  mat A   = join_vert(B, join_horiz(eye<mat>(c-r,c-r), zeros<mat>(c-r,r)));
  return(A);
}

// [[Rcpp::export]]
Rcpp::Date last_year_holiday(Rcpp::Date today,
                             Rcpp::Date this_holiday,
                             Rcpp::Date last_holiday){
  int ddiff = today-this_holiday;
  Date h_date = last_holiday + ddiff;
  return h_date; // 2000-01-02
}

// [[Rcpp::export]]
Rcpp::Date last_year(Rcpp::Date today){
  Date d(today.getYear()-1,today.getMonth() ,today.getDay());
  //int ddiff = today-d;
  return d; // 2000-01-02
}

// [[Rcpp::export]]
Rcpp::Date years_ago(Rcpp::Date today,
                     int shift){
  Date d(today.getYear()-shift,today.getMonth() ,today.getDay());
  //int ddiff = today-d;
  return d; // 2000-01-02
}

// [[Rcpp::export]]
Rcpp::Date next_year_holiday(Rcpp::Date today,
                             Rcpp::Date this_holiday,
                             Rcpp::Date next_holiday){
  int ddiff = today-this_holiday;
  Date h_date = next_holiday + ddiff;
  return h_date; // 2000-01-02
}

// [[Rcpp::export]]
Rcpp::Date next_year(Rcpp::Date today){
  Date d(today.getYear()+1,today.getMonth() ,today.getDay());
  //int ddiff = today-d;
  return d; // 2000-01-02
}

// [[Rcpp::export]]
arma::uword which_date_leq(Rcpp::Date date,
                           std::vector<Date> Dvec){
  uword j;
  for(j = 0; j < Dvec.size(); j++){
    if(Dvec[j] > date){
      break;
    }
  }
  return(j); //R indexing
}

// [[Rcpp::export]]
arma::uword which_date_geq(Rcpp::Date date,
                           std::vector<Date> Dvec){
  uword j;
  for(j = Dvec.size(); j>0; j--){
    if(Dvec[j-1] < date){
      break;
    }
  }
  return(j+1); //R indexing
}

// [[Rcpp::export]]
arma::uword which_date_closest(Rcpp::Date date,
                               std::vector<Date> Dvec){
  uword j;
  for(j = 0; j < Dvec.size(); j++){
    if(Dvec[j] > date){
      break;
    }
  }
  if(j > 1 && j < Dvec.size()){
    int d_less = date - Dvec[j-1];
    int d_more = Dvec[j] - date;
    if(d_less > d_more){
      j++;
    }
  }
  return(j); //R indexing
}

// [[Rcpp::export]]
arma::uvec which_date_closest_ordered(std::vector<Date> FromVec,    //for each date in this vector
                                      std::vector<Date> IndVec){    //find the index of the closest date in this vector
  int d_more, d_less;
  uword j = 0;
  uword k = 1;
  uvec idx(FromVec.size(), fill::ones);
  idx = idx*IndVec.size();
  while(j<FromVec.size()){
    while(k<IndVec.size()){
      if(FromVec[j] < IndVec[k]){
        d_less = FromVec[j] - IndVec[k-1];
        d_more = IndVec[k] - FromVec[j];
        if(d_more<d_less){
          idx(j) = k + 1; //R indexing
          k++;
          break;
        }else{
          idx(j) = k; //R indexing
          break;
        }
      }else{
        k++;
      }
    }
    j++;
  }
  return(idx);
}

// [[Rcpp::export]]
arma::uword day(Rcpp::Date date){
  uword out = date.getDay();
  return(out);
}

// [[Rcpp::export]]
Rcpp::Date replace_day(Rcpp::Date date,
                        int new_day){
  Date d(date.getYear(),date.getMonth(),new_day);
  return(d);
}


// [[Rcpp::export]]
int MonthDays(double year,
              double month){
  int days;
  if((month == 1) || (month == 3) || (month == 5) || (month == 7) || (month == 8) || (month == 10) || (month == 12)){
    days = 31;
  }
  else if((month == 4) || (month == 6) || (month == 9) || (month == 11) ){
    days = 30;
  }
  else if(round((year-1940)/4) == ((year-1940)/4) ){
    days = 29;
  }
  else{
    days = 28;
  }
  return(days);
}

//return last day for the given month
// [[Rcpp::export]]
std::vector<Date> end_of_month(std::vector<Date> date){
  std::vector<Date> d(date.size());
  Rcpp::Date tmp;
  int days;
  for(uword j=0; j<date.size(); j++){
    tmp  = date[j];
    days = MonthDays(tmp.getYear(), tmp.getMonth());
    d[j] = Date(tmp.getYear(), tmp.getMonth(), days);
  }
  return(d);
}

//return last day for the given month
// [[Rcpp::export]]
Rcpp::Date end_of_month_date(Rcpp::Date date){
  int days = MonthDays(date.getYear(), date.getMonth());
  Rcpp::Date d(date.getYear(), date.getMonth(), days);
  return(d);
}

//return first day for the given month
// [[Rcpp::export]]
std::vector<Date> first_of_month(std::vector<Date> date){
  std::vector<Date> d(date.size());
  Rcpp::Date tmp;
  for(uword j=0; j<date.size(); j++){
    tmp  = date[j];
    d[j] = Date(tmp.getYear(), tmp.getMonth(), 1);
  }
  return(d);
}

// [[Rcpp::export]]
std::vector<Date> pseudo_weekly_sequence(Rcpp::Date start,
                                        arma::uword length){
  uword year = start.getYear();
  uword month = start.getMonth();
  uword day = start.getDay();
  uword d = floor(day/7) - 1;
  uword zip = 0;
  d = std::max(zip,d);
  uvec days;
  days << 7 << 14 << 21 << 28 << endr;
  std::vector<Date> out(length);
  uword j = 0;
  while(j<length){
    month = 1;
    while(month<13){
      d = 0;
      while(d<4){
        if(d<3){
          out[j] = Date(year,month,days[d]);
        }else{
          out[j] = end_of_month_date(Date(year,month,days[d]));
        }
        d++;
        j++;
      }
      month++;
    }
    year++;
  }
  return(out);
}
//std::vector<Date>
// [[Rcpp::export]]
Rcpp::Date pseudo_weekly_date(Rcpp::Date date){
  double tmp = date.getDay();
  int week = std::ceil(tmp/7);
  int day;
  if(week==4 || week == 5){
    day = MonthDays(date.getYear(), date.getMonth());
  }else{
    day = week*7;
  }
  Rcpp::Date d(date.getYear(), date.getMonth(), day);
  return(d);
}

// std::vector<Date> d(date.size());
// Rcpp::Date tmp;
// int day;
// double td;
// int week;
// for(uword j=0; j<date.size(); j++){
//   tmp  = date[j];
//   td = tmp.getDay()/7;
//   week = std::ceil(td);
//   if(week==4 || week == 5){
//     day = MonthDays(tmp.getYear(), tmp.getMonth());
//   }else{
//     day = week*7;
//   }
//   d[j] = Date(tmp.getYear(), tmp.getMonth(), day);


// [[Rcpp::export]]
List SARMA(arma::vec Y, //univariate data
               arma::mat p, //AR parameters
               arma::mat q, //MA parameters
               arma::mat P, //seasonal AR parameters
               arma::mat Q, //seasonal MA parameters
               arma::umat P_lag, //seasonal AR lags
               arma::umat Q_lag){ //seasonal MA lags


  uword T   = Y.n_rows;
  uword T_long = std::max(Q_lag.n_rows,P_lag.n_rows);
  T_long = std::max(T,T_long);
  uword sp  = p.n_cols;
  uword sq  = q.n_cols;
  uword pq  = std::max(sp,sq);
  double s_AR, s_MA;
  vec E(T, fill::zeros);
  vec eps(T, fill::zeros);
  vec seas(T_long, fill::zeros);
  vec YP(T, fill::zeros);
  vec Y_sa(T, fill::zeros);
  uvec all_ind;
  //where to start the iterations
  uvec Pl(2, fill::zeros);
  uvec Ql(2, fill::zeros);
  if(P.n_cols>0){
    Pl = ind2sub(size(P_lag), max(find(P_lag==0)));
  }
  if(Q.n_cols>0){
    Ql = ind2sub(size(Q_lag), max(find(Q_lag==0)));
  }
  uword srt = std::max(pq + Pl[0], pq + Ql[0]);

  if(sp==0 && sq==0){
    throw("Seasonal adjustment requires at least 1 AR lag");
  }else if(sp>0 && sq==0){
    if(P.n_cols==0 && Q.n_cols==0){
      for(uword t = srt; t<T; t++){
        if( !Y(span(t-sp,t)).is_finite() ) continue;
        YP(t) = as_scalar(p*flipud(Y(span(t-sp,t-1)))); //adding non-seasonal components
        E(t) = Y(t) - YP(t);
      }
    }else if(P.n_cols>0 && Q.n_cols==0){
      for(uword t = srt; t<T; t++){
        s_AR = as_scalar(P*Y(P_lag.row(t))); //seasonal AR covariates
        if( !Y(span(t-sp,t)).is_finite() || !is_finite(s_AR) ) continue;
        seas(t) = s_AR; //seasonal component
        Y_sa(t) = Y(t) - seas(t);
        YP(t) = as_scalar(p*flipud(Y_sa(span(t-sp,t-1)))); //adding non-seasonal components
        E(t) = Y_sa(t) - YP(t);
      }
    }else if(P.n_cols==0 && Q.n_cols>0){
      for(uword t = srt; t<T; t++){
        s_MA = as_scalar(Q*E(Q_lag.row(t))); //seasonal MA covariates
        if( !Y(span(t-sp,t)).is_finite() || !is_finite(s_MA) ) continue;
        seas(t) = s_MA; //seasonal component
        Y_sa(t) = Y(t) - seas(t);
        YP(t) = as_scalar(p*flipud(Y_sa(span(t-sp,t-1)))); //adding non-seasonal components
        E(t) = Y_sa(t) - YP(t);
      }
    }else{
      for(uword t = srt; t<T; t++){
        s_MA = as_scalar(Q*E(Q_lag.row(t))); //seasonal MA covariates
        s_AR = as_scalar(P*Y(P_lag.row(t))); //seasonal AR covariates
        if( !Y(span(t-sp,t)).is_finite() || !is_finite(s_MA) || !is_finite(s_AR) ) continue;
        seas(t) = s_AR + s_MA; //seasonal component
        Y_sa(t) = Y(t) - seas(t);
        YP(t) = as_scalar(p*flipud(Y_sa(span(t-sp,t-1)))); //adding non-seasonal components
        E(t) = Y_sa(t) - YP(t);
      }
    }
  }else if(sp==0 && sq>0){
    throw("Seasonal adjustment requires at least 1 AR lag");
  }else{
    if(P.n_cols==0 && Q.n_cols==0){ //q not used
      throw("Seasonal adjustment with q>0 requires Q>0");
    }else if(P.n_cols>0 && Q.n_cols==0){ //q not used
      throw("Seasonal adjustment with q>0 requires Q>0");
    }else if(P.n_cols==0 && Q.n_cols>0){
      for(uword t = srt; t<T; t++){
        s_MA = as_scalar(Q*eps(Q_lag.row(t))); //seasonal MA covariates
        if( !Y(span(t-pq,t)).is_finite() || !is_finite(s_MA) ) continue;
        seas(t) = s_MA; //seasonal component
        Y_sa(t) = Y(t) - seas(t);
        YP(t) = as_scalar(p*flipud(Y_sa(span(t-sp,t-1)))); //adding MA components
        E(t) = Y_sa(t) - YP(t);
        eps(t) = as_scalar(q*flipud(eps(span(t-sq,t-1)))) + E(t);
      }
    }else{
      for(uword t = srt; t<T; t++){
        s_MA = as_scalar(Q*eps(Q_lag.row(t))); //seasonal MA covariates
        s_AR = as_scalar(P*Y(P_lag.row(t))); //seasonal AR covariates
        if( !Y(span(t-pq,t)).is_finite() || !is_finite(s_MA) || !is_finite(s_AR) ) continue;
        seas(t) = s_AR + s_MA; //seasonal component
        Y_sa(t) = Y(t) - seas(t);
        YP(t) = as_scalar(p*flipud(Y_sa(span(t-sp,t-1)))); //predicting next period SA values
        E(t) = Y_sa(t) - YP(t);
        eps(t) = as_scalar(q*flipud(eps(span(t-sq,t-1)))) + E(t);
      }
    }
  }

  //Get future seasonal adjustments if required
  if(T_long > T){
    if(P.n_cols==0 && Q.n_cols>0){
      for(uword t = T; t<T_long; t++){
        s_MA = as_scalar(Q*eps(Q_lag.row(t))); //seasonal MA covariates
        if( !is_finite(s_MA) ) continue;
        seas(t) = s_MA; //seasonal component
      }
    }else if(P.n_cols>0 && Q.n_cols==0){
      for(uword t = T; t<T_long; t++){
        s_AR = as_scalar(P*Y(P_lag.row(t))); //seasonal AR covariates
        if( !is_finite(s_AR) ) continue;
        seas(t) = s_AR; //seasonal component
      }
    }else{
      for(uword t = T; t<T_long; t++){
        s_MA = as_scalar(Q*eps(Q_lag.row(t))); //seasonal MA covariates
        s_AR = as_scalar(P*Y(P_lag.row(t))); //seasonal AR covariates
        if( !is_finite(s_MA) || !is_finite(s_AR) ) continue;
        seas(t) = s_AR + s_MA; //seasonal component
      }
    }
  } // if(T_long > T)

  //double MSE = as_scalar(trans(E(span(srt,T-1)))*E(span(srt,T-1)))/T;

  uvec ind = find(E); //find non-zero elements of E
  double MSE = as_scalar(trans(E(ind))*E(ind))/(ind.n_elem);

  List Out;
  Out["MSE"] = MSE;
  Out["seas"]  = seas;
  Out["Y_sa"] = Y_sa;
  Out["E"] = E;
  Out["YP"] = YP;
  Out["p"]  = p;
  Out["q"] = q;
  Out["P"] = P;
  Out["Q"] = Q;
  Out["srt"] = srt;

  return(Out);
}

