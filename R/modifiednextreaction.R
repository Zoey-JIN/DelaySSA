#' @title Modified Next Reaction Method Algorithm
#' @description A Stochastic Simulation Algorithm (SSA) without Delays Using Modified Next Reaction Method
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

simulate_reaction_modifiednextreaction <- function(tmax, n_initial, t_initial, S_matrix, k, fun_fr) {
  n_values <- matrix(n_initial)
  t_values <- c(t_initial)
  n <- n_initial
  t <- t_initial
  p_vec <- rep(0, ncol(S_matrix))
  t_vec <- rep(0, ncol(S_matrix))

  f_r <- fun_fr(k,n)
  u2 <- runif(ncol(S_matrix))
  p_vec <- log(1/u2)

  while(t < tmax){
    tau_vec <- (p_vec-t_vec) / f_r
    r <- which.min(tau_vec)
    tau <- tau_vec[r]
    n <- n + S_matrix[,r]

    t_vec <- t_vec+f_r*tau
    u1 <- runif(1)
    p_vec[r] <- p_vec[r]+log(1/u1)
    f_r <- fun_fr(k,n)

    t <- t+tau
    t_values <- c(t_values, t)
    n_values <- cbind(n_values, n)
  }
  return(list(t_values = t_values, n_values = n_values))
}
