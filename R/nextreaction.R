#' @title Next Reaction Method Algorithm
#' @description A Stochastic Simulation Algorithm (SSA) without Delays Using Next Reaction Method
#'
#'
#' @param tmax cutoff time
#' @param n_initial initial species number
#' @param t_initial initial time
#' @param S_matrix the stoichiometric matrix at the initiation time
#' @param k a reaction rate vector
#' @param fun_fr a function representing the propensity function
#'
#' @return the amount of a species and the corresponding time
#' @export
#'
simulate_reaction_nextreaction <- function(tmax, n_initial, t_initial, S_matrix, k, fun_fr) {
  n_values <- matrix(n_initial)
  t_values <- c(t_initial)
  n <- n_initial
  t <- t_initial

  f_r <- fun_fr(k,n)

  u2 <- runif(ncol(S_matrix))
  tau_vec <- -log(u2) / f_r

  while(t < tmax){
    r <- which.min(tau_vec)
    tau <- tau_vec[r]
    n <- n + S_matrix[,r]
    f_r_update <- fun_fr(k,n)
    if(length(tau_vec[-r])>0){
      tau_vec[-r] <- f_r[-r]/f_r_update[-r]*(tau_vec[-r]-tau)+tau
    }
    if(any(is.infinite(tau_vec)|is.na(tau_vec))){
      Inf_index <- which(is.na(is.infinite(tau_vec)|tau_vec))
      tau_vec[Inf_index] <- -log(runif(length(Inf_index))) / f_r_update[Inf_index] + tau
    }
    u1 <- runif(1)
    tau_vec[r] <- 1/f_r_update[r]*log(1/u1)+tau
    f_r <- f_r_update
    t <- tau
    t_values <- c(t_values, t)
    n_values <- cbind(n_values, n)
  }
  return(list(t_values = t_values, n_values = n_values))
}
