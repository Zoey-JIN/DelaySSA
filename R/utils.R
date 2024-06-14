#' @title Delay Time
#'
#' @param element an element in delaytime_list which could be a function or a fixed number
#'
#' @return delaytime which is a fixed number or a stochastic value governed by a function
tau_element <- function(element) {
  if (is.function(element)) {
    return(element())
  } else {
    return(element)
  }
}

#' Delay Effect Matrix
#'
#' @param delay_type vector of reaction type taking on the values 0, 1, or 2
#' @param delaytime_list a list representing the delay time of each reaction
#' @param S_matrix the stoichiometric matrix at the initiation time
#' @param S_matrix_delay the stoichiometric matrix at the completion time
#'
#' @return A 2-row delay_effect_matrix
#' @export
check_delay_relative <- function(delay_type, delaytime_list, S_matrix, S_matrix_delay) {
  zero_check <- apply(S_matrix_delay[, which(delay_type == 0), drop = FALSE], 2, function(col) all(col == 0))

  if (!all(zero_check)) {
    warning("Not all corresponding columns in s_matrix_delay are zeros for indices where delay_type is 0")
  }
  zero_check <- sapply(delaytime_list[which(delay_type == 0)], function(x) x == 0)
  if (!all(zero_check)) {
    warning("Not all corresponding elements in delaytime_list are zeros for indices where delay_type is 0")
  }
  zero_check <- apply(S_matrix[, which(delay_type == 1), drop = FALSE], 2, function(col) all(col == 0))
  if (!all(zero_check)) {
    warning("Not all corresponding elements in s_matrix are zeros for indices where delay_type is 1")
  }

  non_zero_indices <- which(delay_type != 0)
  result_matrix <- matrix(NA, nrow = 2, ncol = 0)

  for (idx in non_zero_indices) {
    match_index <- which(apply(S_matrix, 2, function(col) all(S_matrix_delay[, idx] == col)))
    if (length(match_index) > 0) {
      result_matrix <- cbind(result_matrix, c(match_index[1], idx))
    }
  }
  return(result_matrix)
}

#' @title For Propensity Calculation
#'
#' @param n species number
#' @param reactant_matrx species reactant matrix
#'
#' @return a vector used for calculating propensity function
propensity_n <- function(n,reactant_matrx){
  result <- sapply(1:ncol(reactant_matrx), function(j) {
    prod(sapply(1:nrow(reactant_matrx), function(i) {
      ifelse(reactant_matrx[i,j]==0, 1, Reduce(`*`, (n[i] - (0:(reactant_matrx[i,j]-1)))))
    }))})
  return(result)
}

#' @title The Mean Value of Species i at Time t
#'
#' @param result result of solve_DelaySSA
#' @param i species i
#' @param t time t
#'
#' @return mean value
#' @export
plot_mean <- function(result,i,t){
  n <- lapply(result,function(x) picksample(x,i,t))
  n <- unlist(n)
  return(mean(n))
}

#' @title Probability Density Functions Conversion
#'
#' @param a_vector a vector
#'
#' @return PDF
#' @export
convert_pdf <- function(a_vector) {
  freq <- table(a_vector)
  prob <- freq / sum(freq)
  unique_values <- as.integer(names(prob))
  output <- list(unique_values, prob)
  return(output)
}

#' @title Results Sampling
#'
#' @param list_output a outcome of a single simulation
#' @param i species i
#' @param t time t
#'
#' @return species number
#' @export
picksample <- function(list_output,i=1,t){
  index <- which.min(list_output$t_values<t)-1
  n_output <- list_output$n_values[i,index]
  return(n_output)
}


#' @title A Stochastic Simulation Algorithm (SSA) with Delays
#' @description DelaySSA implements a stochastic simulation algorithm (SSA) with delays in R through different algorithms.
#'
#'
#' @param algorithm  "DelayMNR" (default), "DelayReject", or "DelayDirect". If without delay, "Direct", "MNR" and "NR"  are recommended and the parameters of "S_matrix_dalay" "delay_type" and "delaytime_list" can be omitted.
#' @param sample_size the number of repetitions
#' @param tmax cutoff time
#' @param n_initial initial species number
#' @param t_initial initial time
#' @param S_matrix the stoichiometric matrix at the initiation time
#' @param S_matrix_delay the stoichiometric matrix at the completion time
#' @param k a reaction rate vector
#' @param reactant_matrx species reactant matrix
#' @param delay_type the reaction type vector taking on the values 0, 1, or 2
#' @param delaytime_list a list representing the delay time of each reaction
#'
#' @return A list contains sublists including the amount of a species and the corresponding time
#' @export
simulation_DelaySSA <- function(algorithm = "DelayMNR", sample_size, tmax, n_initial, t_initial, S_matrix, S_matrix_delay = NULL, k, reactant_matrx, delay_type = NULL , delaytime_list = NULL) {
    algorithm_chosen <- algorithm
    sample <- sample_size
    if (!(is.null(delay_type) || is.null(delaytime_list) || is.null(S_matrix_delay))) {
      delay_effect_matrix <- check_delay_relative(delay_type, delaytime_list, S_matrix, S_matrix_delay)}
    fun_fr <- function(k,n){
      if (is.function(k)) {
        k_mask <- k(n)
      }else {
        k_mask <- k
      }
      return(k_mask*propensity_n(n,reactant_matrx))
    }
    result <- switch(algorithm_chosen,
                    "DelayDirect" = lapply(1:sample, function(x) simulate_reaction_delay_direct(tmax=, n_initial, t_initial, S_matrix, S_matrix_delay, k,fun_fr, delay_type, delaytime_list, delay_effect_matrix)),
                    "DelayMNR" = lapply(1:sample, function(x) simulate_reaction_delay_modifiednextreaction(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix)),
                    "DelayRejection" = lapply(1:sample, function(x) simulate_reaction_delay_rejection(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix)),
                    "Direct" = lapply(1:sample, function(x) simulate_reaction(tmax, n_initial, t_initial, S_matrix, k, fun_fr)),
                    "MNR" = lapply(1:sample, function(x) simulate_reaction_modifiednextreaction(tmax, n_initial, t_initial, S_matrix, k, fun_fr)),
                    "NR" = lapply(1:sample, function(x) simulate_reaction_nextreaction(tmax, n_initial, t_initial, S_matrix, k, fun_fr)),
                    "Error: No such algorithm!"
                    )
   return(result)
}

#' @title Density at Time t_pick
#'
#' @param result result of solve_DelaySSA
#' @param t_pick time t_pick
#'
#' @return a picture showing the density at t_pick
#' @export
plot_SSA_density <- function(result,t_pick){
  num_columns <- nrow(result[[1]]$n_values)
  data_list <- lapply(1:num_columns, function(i) {
    n <- lapply(result, function(x) picksample(x, i, t_pick))
    n <- unlist(n)
    plot_xy <- convert_pdf(n)
    if(!all(data.frame(percentage = plot_xy)[,1]== data.frame(percentage = plot_xy)[,2])){
      warning("Error in Calculating Density Table")
    }
    data.frame(quantity = data.frame(percentage = plot_xy)[,1], percentage = data.frame(percentage = plot_xy)[,3], Specie = paste("Specie", i))
  })
  plot_data <- do.call(rbind, data_list)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = quantity, y = percentage, color = Specie)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Quantity", y = "Percentage", title = paste("PDF of Specie at time",t_pick)) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_blank()
    )

}

#' @title Mean of Species over Time
#'
#' @param result result of solve_DelaySSA
#' @param t time series
#' @param n_initial initial species number
#' @param t_initial initial time
#'
#' @return a picture showing the mean through time
#' @export
plot_SSA_mean <- function(result,t = seq(0, tmax, by = 1),n_initial=n_initial,t_initial=t_initial){
  num_columns <- nrow(result[[1]]$n_values)
  data_list <- lapply(1:num_columns, function(i) {
    n <- lapply(t,function(x) plot_mean(result,i,x))
    n <- unlist(n)
    if(t[1]==t_initial)
      n[1] <- n_initial[i,]
    data.frame(t = t, quantity = n, Specie = paste("Specie", i))
  })
  plot_data <- do.call(rbind, data_list)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = t, y = quantity, color = Specie)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "T", y = "Quantity", title = "Mean of Specie") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_blank()
    )
}
#' @title Propensity Function
#'
#' @param k reaction rate
#' @param n species number
#'
#' @return Propensity Function
fun_fr <- function(k,n){
  if (is.function(k)) {
    k_mask <- k(n)
  }else {
    k_mask <- k
  }
  return(k_mask*propensity_n(n,product_list))
}
