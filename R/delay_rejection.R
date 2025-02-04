#' @title Delay Rejection Method Algorithm
#' @description A Stochastic Simulation Algorithm (SSA) with Delays Using Rejection Method
#'
#'
#' @param tmax cutoff time
#' @param n_initial initial species number
#' @param t_initial initial time
#' @param S_matrix the stoichiometric matrix at the initiation time
#' @param S_matrix_delay the stoichiometric matrix at the completion time
#' @param k a reaction rate vector
#' @param fun_fr a function representing the propensity function
#' @param delay_type the reaction type vector taking on the values 0, 1, or 2
#' @param delaytime_list a list representing the delay time of each reaction
#' @param delay_effect_matrix a matrix representing that reaction without delay affects reaction with delay
#' @param reactant_matrix_delay species reactant matrix in delay part
#'
#' @return the amount of a species and the corresponding time
#' @export
#'

simulate_reaction_delay_rejection <- function(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix, reactant_matrix_delay) {
  n_values <- matrix(n_initial)
  t_values <- c(t_initial)
  n <- n_initial
  t <- t_initial

  Tstruct <- vector("list", length = 2)
  # # The first line represents the time, and the second line indicates the reaction.
  while (t < tmax) {
    # tau
    u1 <- runif(1)
    f_r <- fun_fr(k,n)
    lambda_sum <- sum(f_r)
    tau <- -log(u1) / lambda_sum

    if(any(Tstruct[[1]]<t+tau)){
      t <- Tstruct[[1]][1]
      r <- Tstruct[[2]][1]
      if(delay_type[r]==0){
        print("warning")
        # n <- n + S_matrix[,r]
      } else if (delay_type[r]==1) {
        n <- n + S_matrix_delay[,r] + reactant_matrix_delay[,r]
      } else if (delay_type[r]==2) {
        n <- n + S_matrix_delay[,r]
      } else {
        print("Stop: wrong with the reaction type")
      }
      Tstruct[[1]] <- Tstruct[[1]][-1]
      Tstruct[[2]] <- Tstruct[[2]][-1]
    }else{
      u2 <- runif(1)
      r <- which.max(cumsum(f_r) > u2 * lambda_sum)
      t <- t+tau
      if(delay_type[r]==0){
        if(r %in% delay_effect_matrix[1,]){
          effect_r <- delay_effect_matrix[2,][which(delay_effect_matrix[1, ] == r)]
          for (element in effect_r) {
            drop_index <- sample(which(Tstruct[[2]] == element), 1) 
            if(delay_type[Tstruct[[2]][drop_index]]==1){
              n <- n + reactant_matrix_delay[,Tstruct[[2]][drop_index]]
            }
            Tstruct[[1]] <- Tstruct[[1]][-drop_index]
            Tstruct[[2]] <- Tstruct[[2]][-drop_index]
          }
        }
        n <- n + S_matrix[,r]
      } else if (delay_type[r]==1) {
        n <- n + S_matrix[,r]
        n <- n - reactant_matrix_delay[,r]
        add_tau <- tau_element(delaytime_list[[r]])+t
        index <- findInterval(add_tau,Tstruct[[1]])
        if (index == 0) {
          Tstruct[[1]] <- c(add_tau, Tstruct[[1]])
          Tstruct[[2]] <- c(r, Tstruct[[2]])
        } else if (index == length(Tstruct[[1]])) {
          Tstruct[[1]] <- c(Tstruct[[1]], add_tau)
          Tstruct[[2]] <- c(Tstruct[[2]], r)
        } else {
          Tstruct[[1]] <- c(Tstruct[[1]][1:index], add_tau, Tstruct[[1]][(index + 1):length(Tstruct[[1]])])
          Tstruct[[2]] <- c(Tstruct[[2]][1:index], r, Tstruct[[2]][(index + 1):length(Tstruct[[2]])])
        }
      } else if (delay_type[r]==2) {
        n <- n + S_matrix[,r]
        add_tau <- tau_element(delaytime_list[[r]])+t
        index <- findInterval(add_tau,Tstruct[[1]])
        if (index == 0) {
          Tstruct[[1]] <- c(add_tau, Tstruct[[1]])
          Tstruct[[2]] <- c(r, Tstruct[[2]])
        } else if (index == length(Tstruct[[1]])) {
          Tstruct[[1]] <- c(Tstruct[[1]], add_tau)
          Tstruct[[2]] <- c(Tstruct[[2]], r)
        } else {
          Tstruct[[1]] <- c(Tstruct[[1]][1:index], add_tau, Tstruct[[1]][(index + 1):length(Tstruct[[1]])])
          Tstruct[[2]] <- c(Tstruct[[2]][1:index], r, Tstruct[[2]][(index + 1):length(Tstruct[[2]])])
        }
      } else {
        print("Stop: wrong with the reaction type")
      }
    }
    if(t<tail(t_values,1)){
      break
    }
    t_values <- c(t_values, t)
    n_values <- cbind(n_values, n)
  }
  return(list(t_values = t_values, n_values = n_values))
}
