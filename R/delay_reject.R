#' @title Delay Rejection Method Algorithm
#' @description A Stochastic Simulation Algorithm (SSA) with Delays Using Rejection Method
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

simulate_reaction_delay_rejection <- function(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, fun_fr, delay_type, delaytime_list, delay_effect_matrix) {
  n_values <- matrix(n_initial)
  t_values <- c(t_initial)
  n <- n_initial
  t <- t_initial

  Tstruct <- vector("list", length = 2)
  # 第一行为时间 第二行为第几个反应
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
        n <- n + S_matrix[,r]
      } else if (delay_type[r]==1) {
        n <- n + S_matrix_delay[,r]
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
    }
    t_values <- c(t_values, t)
    n_values <- cbind(n_values, n)
  }
  return(list(t_values = t_values, n_values = n_values))
}
