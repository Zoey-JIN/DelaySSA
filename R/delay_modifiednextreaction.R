#' @title Delay Modified Next Reaction Method Algorithm
#' @description A Stochastic Simulation Algorithm (SSA) with Delays Using Modified Next Reaction Method
#'
#'
#' @param tmax cutoff time
#' @param n_initial initial species number
#' @param t_initial initial time
#' @param S_matrix the stoichiometric matrix at the initiation time
#' @param S_matrix_delay the stoichiometric matrix at the completion time
#' @param fun_fr a function representing the propensity function
#' @param delay_type vector of reaction type taking on the values 0, 1, or 2
#' @param delaytime_list a list representing the delay time of each reaction
#' @param delay_effect_matrix a matrix representing that reaction without delay affects reaction with delay
#'
#' @return the amount of a species and the corresponding time
#' @export
#'

simulate_reaction_delay_modifiednextreaction <- function(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, fun_fr, delay_type, delaytime_list, delay_effect_matrix) {
  n_values <- matrix(n_initial)
  t_values <- c(t_initial)
  n <- n_initial
  t <- t_initial

  p_vec <- rep(0, ncol(S_matrix))
  t_vec <- rep(0, ncol(S_matrix))
  Tstruct <- vector("list", length = 2)

  f_r <- fun_fr(k,n)
  u2 <-  runif(ncol(S_matrix))
  p_vec <- log(1/u2)

  while (t < tmax) {

    tau_vec <- (p_vec-t_vec)/f_r + t

    min_1 <- min(tau_vec)
    r_1 <- which.min(tau_vec)
    min_2 <- Tstruct[[1]][1]
    r_2 <-  Tstruct[[2]][1]
    if(is.na(min_2) || length(min_2) == 0 || min_1 < min_2){
      r <- r_1
      tau <- min_1-t
      t <- min_1
      if(delay_type[r]==0){
        if(r %in% delay_effect_matrix[1,]){
          effect_r <- delay_effect_matrix[2,][which(delay_effect_matrix[1, ] == r)]
          drop_index <- sample(which(Tstruct[[2]] == effect_r),1)
          Tstruct[[1]] <- Tstruct[[1]][-drop_index]
          Tstruct[[2]] <- Tstruct[[2]][-drop_index]
        }
        n <- n + S_matrix[,r]
      } else if (delay_type[r]==1) {
        add_tau <- tau_element(delaytime_list[[r]])+t
        index <- findInterval(add_tau,Tstruct[[1]])
        Tstruct[[1]] <- append(Tstruct[[1]], add_tau, after = index)
        Tstruct[[2]] <- append(Tstruct[[2]], r, after = index)
      } else if (delay_type[r]==2) {
        n <- n + S_matrix[,r]
        add_tau <- tau_element(delaytime_list[[r]])+t
        index <- findInterval(add_tau,Tstruct[[1]])
        Tstruct[[1]] <- append(Tstruct[[1]], add_tau, after = index)
        Tstruct[[2]] <- append(Tstruct[[2]], r, after = index)
      } else {
        print("Stop: wrong with the reaction type")
      }
      u1 <- runif(1)
      p_vec[r] <- p_vec[r]+log(1/u1)
    } else{
      r <- r_2
      tau <- min_2-t
      t <- min_2
      n <- n+S_matrix_delay[,r]
      Tstruct[[1]] <- Tstruct[[1]][-1]
      Tstruct[[2]] <- Tstruct[[2]][-1]
    }

    t_vec<- t_vec+f_r*tau

    f_r <- fun_fr(k,n)
    t_values <- c(t_values, t)
    n_values <- cbind(n_values, n)
  }

  return(list(t_values = t_values, n_values = n_values))
}
