
#' Get B-spline basis function with a fixed gap between knots
#' 
#' @param xx Covariate for which spline basis is to be defined
#' @param knot_gap The space between adjacent knots, using the same units as 
#' \code{xx}
#' 
#' @return A matrix of basis functions, with three knots outside the range of the 
#' data at each end.
get_static_basis_function <- function(xx, knot_gap, deriv=F){
  max_x <- max(xx)
  first_x <- min(xx)
  knot_locations <-  seq(first_x - (knot_gap*3), 
                         max_x + 1 + knot_gap*3, 
                         knot_gap)
  basis <- splines::splineDesign(knot_locations, xx, outer.ok=T, derivs = deriv)
  return(basis)
}



#' Return a penalty matrix
#' 
#' @param n_basis the number of bases to be penalised.
#' @param age_penalty_order The order of the penalty
#' 
#' @return A matrix that can be used to penalise differences.
make_penalty_matrix <- function(n_basis, age_penalty_order){
  dd <- diff(diag(n_basis), differences = age_penalty_order)
  penalty_matrix <- t(dd) %*% dd
  return(penalty_matrix)
}


#' Find a penalty for the null space of an existing penalty matrix. See
#' \code{mgcv} package and the documentation for jagam.
#' 
#' @param penalty_matrix An existing penalty matrix, for example from 
#' \code{make_penalty_matrix}
#' 
#' @return A matrix that penalises the null space of \code{penalty_matrix}
make_null_penalty <- function(penalty_matrix){
  pen_eig <- eigen(penalty_matrix)
  inds <- which(abs(pen_eig$values)<1e-10)
  u0 <- pen_eig$vector[,inds]
  penalty_matrix_null <- u0 %*% t(u0)
  return(penalty_matrix_null)
}


#' Get the conditional covariance matrix for the period parameters,
#'
#' A distribution for the period innovations is needed for model fitting, 
#' conditional on the identifiability constraints. 
#' This function constructs the covariance matrix of this 
#' distribution for use in model fitting. This covariance matrix may be scaled 
#' by a variance parameter during the sampling process.
#'
#' @param n_years The number of years in the data.
#' @param n_forecast_years Number of forecast years to include in the covariance
#' matrix. Only useful in the two-sex case.
#' @param growth_constraint Is the period effect constrained to show no linear
#' growth
#' @param quad_constraint Is the period effect constrained to show no quadratic
#' growth

#' 
#' @return The raw covariance matrix for period effect, conditional on the 
#' constraints
get_period_conditional_cov_matrix <- function(n_years, 
                                              n_forecast_years=0,
                                              growth_constraint=T,
                                              quad_constraint=F,
                                              basis_mat=diag(n_years + 
                                                             n_forecast_years)){
  
  constraint <- get_constraint(basis_mat,
                               sum_zero = T, 
                               zero_growth=growth_constraint,
                               zero_quad = quad_constraint,
                               n_forecast_years = n_forecast_years)
  # get cumulative sum matrix
  S <- solve(get_difference_matrix(n_years + n_forecast_years))
  CS <- constraint %*% S
  # Which bases should be constrained 
  n_constraint <- dim(constraint)[1]
  con_ind <- 1:n_constraint
  
  result <- get_conditional_cov_matrix(CS, con_ind)
  return(result)
}

#' Construct a conditional covariance matrix
#' 
#' Given an unconditional covariance matrix, and the indexes to which the 
#' constraints apply, construct an unscaled conditional covariance matrix, based
#' on standard multivariate normal results. 
#' 
#' @param CS Unconditional covariance matrix
#' @param con_ind The indexes to which the constraints apply.
#'
get_conditional_cov_matrix <- function(CS, con_ind){
  ZZ <- get_constraint_transformation_matrix(CS, con_ind)
  Sigma <- ZZ %*% t(ZZ)
  Tau <- (Sigma[-con_ind,-con_ind] - Sigma[-con_ind,con_ind] %*%
            solve(Sigma[con_ind,con_ind]) %*%
            Sigma[con_ind,-con_ind])
  inv_constraint <- solve(ZZ)[,-con_ind]
  return(list(Tau=Tau, inv_constraint=inv_constraint))
}


#' Construct matrix to allow conditioning on constraints holding.
#' 
#' Constructs matrix that transforms from the initial differenced parameter 
#' space to one where some rows are zero if the constraints on the cumulative 
#' sums of the parameters (specifed in `constraint`) hold, allow conditioning on
#' the constraints holding true.
#' 
#' @param constraint  A matrix where each row is a constraint which must equal 
#' zero.
#' 
#' @return A square matrix for which, for $k$ constraints, will contain $k$  
#' rows that will sum to zero if the desired constraints on the cumulative sums 
#' hold.
#'
get_constraint_transformation_matrix <- function(constraint, con_ind){
  n_basis <- dim(constraint)[2]
  ZZ <- diag(n_basis)
  ZZ[con_ind,] <- constraint
  return(ZZ)
}


#' Construct constraints for a basis matrix
#' 
#' Sometimes the output of a basis function might need to be constrained. This 
#' function constructs the constraint matrix needed to achieves this.
#' Combinations of constraints can be specificed. 
#' 
#' @param basis_matrix The basis matrix to be constrained
#' @param sum_zero Switch indicating whether sum-to-zero constraint should hold
#' @param zero_growth Switch indicating whether zero-growth constraint should 
#' hold
#' @param zero_quad Switch indicating whether zero-quadratic-growth constraint 
#' should hold
#' @param diff_sum_zero Switch indicating whether a sum-to-zero constraint 
#' should hold on the differenced scale
#' @param diff_zero_growth Switch indicating whether a zero-growth constraint 
#' should hold on the differenced scale
#' @param initial_zero Switch indicating whether the first value of the smooth 
#' function should be set equal to zero.
#' @param initial_value_equal Switch indicating whether the first and last 
#' values should be equal
#' @param first_last_equal Switch indicating whether the first and last values
#' should be set to be equal.
#' @param last_zero Switch indicating whether the first and last values
#' should be set to be equal.
#' @param n_forecast how many years of forecasts are there, for which we do not
#' want constraints to hold
#' @param zero_pad number of zero to prepend to the constraint.
#' @return A constraint matrix for use in constructing a constrained basis
#' @export
#'
get_constraint <- function(basis_matrix, sum_zero=F, zero_growth=F, zero_quad=F,
                           diff_sum_zero=F, diff_zero_growth=F, 
                           initial_zero=F, initial_value_equal=F,
                           first_last_equal=F, last_zero=F,
                           n_forecast_years=0, zero_pad = 0){
  cons <- list()
  n_elems <- dim(basis_matrix)[1] - n_forecast_years - zero_pad
  n_basis <- dim(basis_matrix)[2]
  if (sum_zero){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           rep(1,n_elems), 
                           rep(0,n_forecast_years)) 
                         %*% basis_matrix))
  }
  if (zero_growth){
    cons <- c(cons ,list(c(rep(0,zero_pad), seq(- (n_elems - 1)/2,
                                                (n_elems - 1)/2),
                           rep(0,n_forecast_years))  %*% basis_matrix)  )
  }
  if (diff_zero_growth){
    cons <- c(cons, list(c(rep(0,zero_pad), 
                           (n_elems - 2) / 2, rep(-1, n_elems - 2),
                           (n_elems - 2) /2, 
                           rep(0, n_forecast_years)) %*% basis_matrix))
  }
  if (diff_sum_zero){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           -1, rep(0, n_elems - 2), 1, rep(0,n_forecast_years))
                         %*% basis_matrix ))
  }
  if (zero_quad){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           seq(-(n_elems - 1)/2, (n_elems - 1)/2)**2,
                           rep(0, n_forecast_years)) %*% basis_matrix))
  }
  if (initial_zero){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           1, rep(0, n_elems + n_forecast_years - 1)) 
                         %*% basis_matrix))
  }
  if (initial_value_equal){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           1,-1, rep(0, n_basis + n_forecast_years - 2))))
  }
  if (first_last_equal){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           1,rep(0, n_elems - 2),-1, rep(0,n_forecast_years)) 
                         %*% basis_matrix))
  }
  if (last_zero){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           rep(0, n_elems-1), 1,
                           rep(0,n_forecast_years))
                         %*% basis_matrix))
  }
  
  C <- t(sapply(cons, function(x) x))
  return(C)
}


#' Construct the first-difference matrix operator
#' 
#' Returns square matrix $D$ with rows and columns $n$ that computes the first
#' difference of anything it is post-multiplied by.
#' 
#' @param n The number of rows in the difference matrix.
#' 
#' @return An $n$ x $n$ matrix with 1 on the diagonal and -1 on the first 
#' lower off-diagnoal
get_difference_matrix <- function(n){
  D <- matrix(0, n, n)
  D[2:n, 1:(n - 1)] <- - diag(n - 1)
  D <- D + diag(n)
  return(D)
}

# take an arbitrary integer vector and make it start from 1.
make_index <- function(xx){
  xx - min(xx) + 1
}

#' Construct data for input into stan models
#' 
#' @param rate_df Data from the Human Mortality Database containing Deaths and 
#' Exposures for multiple countries.
#' @param config A named list describing properties of the model to be estimated
#'
#' @return A named list containing inputs to the stan model defined by config
make_stan_data <- function(rate_df, config){
  rate_df %<>% filter(Exposure!=0) %>%
    mutate(country_f = as.integer(as.factor(Country)))
   
  N <- dim(rate_df)[1]
  ages <- rate_df$Age %>% unique()
  time_index <- make_index(rate_df$Year)
  age_index <- make_index(rate_df$Age)
  country_index <- make_index(rate_df$country_f)
  country_decode <- rate_df %>% select(Country, country_f) %>% 
    unique()
  n_countries <- length(unique(rate_df$Country))
  
  growth_constraint <- T
  quad_constraint <- F
  if ("quad_constraint" %in% names(config)){
    quad_constraint <- config$quad_constraint
  }
  
  n_period_constraint <- 1 + growth_constraint + quad_constraint 
  # period
  n_years <- length(unique(time_index))
  n_forecast_years <- config$n_forecast_years
  period_cov <- get_period_conditional_cov_matrix(n_years, 
                                            n_forecast_years,
                                            growth_constraint=growth_constraint,
                                            quad_constraint=quad_constraint)
  L_k <- t(chol(period_cov$Tau))
  inv_constraint <- period_cov$inv_constraint
  
  pen_order <- config$pen_order
  age_basis <- get_static_basis_function(ages, config$knot_gap_frontier)
  age_basis_re <- get_static_basis_function(ages, config$knot_gap_countries)
  
  n_basis <- dim(age_basis)[2]
  n_basis_re <- dim(age_basis_re)[2]
  
  penalty_age <- make_penalty_matrix(n_basis, pen_order)
  penalty_age_null  <- make_null_penalty(penalty_age)
  penalty_age_re <- make_penalty_matrix(n_basis_re, pen_order)
  penalty_age_re_null  <- make_null_penalty(penalty_age_re)
  
  
  
  stan_data <- list(n_ages=length(ages),
                    N=N,
                    n_countries = n_countries,
                    n_basis=n_basis,
                    n_years = length(unique(time_index)),
                    n_forecast_years=n_forecast_years,
                    n_basis_re=n_basis_re,
                    n_period_constraint=n_period_constraint,
                    age_basis=age_basis,
                    age_basis_re=age_basis_re,
                    penalty_age=penalty_age,
                    penalty_age_re=penalty_age_re,
                    penalty_age_null=penalty_age_null,
                    penalty_age_re_null=penalty_age_re_null,
                    deaths=round(rate_df$Deaths),
                    expos=rate_df$Exposure,
                    age_index=age_index,
                    country_index=country_index,
                    time_index=time_index,
                    country_decode=country_decode,
                    L_k=L_k,
                    inv_constraint=inv_constraint
                    
  )
  return(stan_data)
}

