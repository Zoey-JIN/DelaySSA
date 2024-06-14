#' @title Gillespie Algorithm
#' @description A Stochastic Simulation Algorithm (SSA) without Delays Using Gillespie Method
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
simulate_reaction <- function(tmax, n_initial, t_initial, S_matrix, k, fun_fr) {
  n_values <- matrix(n_initial)
  t_values <- c(t_initial)
  n <- n_initial
  t <- t_initial

  while (t < tmax) {
    cumulative_sum <- 0
    u1 <- runif(1)
    u2 <- runif(1)

    # f(n)
    f_r <- fun_fr(k,n)
    lambda_sum <- sum(f_r)
    tau <- -log(u1) / lambda_sum

    r <- which.max(cumsum(f_r) > u2 * lambda_sum)

    n <- n + S_matrix[,r]
    t <- t + tau
    t_values <- c(t_values, t)
    n_values <- cbind(n_values, n)
  }

  return(list(t_values = t_values, n_values = n_values))
}
