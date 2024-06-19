# DelaySSA
DelaySSA is an R package designed to simulate non-Markovian models of gene expression using a delayed the Gillespie algorithm. This package is particularly useful for researchers and practitioners in computational biology and systems biology who are interested in modeling the dynamics of gene expression with intrinsic delays.

## Features
- Implementation of the delayed Gillespie algorithm [1-4].
- Support for modeling gene expression with various delay distributions.
- Easy-to-use interface for defining reactions with delays.
- Capability to handle complex gene regulatory networks.
- Visualization tools for analyzing simulation results.

## Installation
Download the installation package from GitHub
```
devtools::install_github("Zoey-JIN/DelaySSA")
```

Or download and install locally
```
devtools::install("~/DelaySSA")
```

Then load the package
```
library("DelaySSA")
```

## Example
Molecular $S_1$ binds $S_2$ and then disappear with the reaction rate $k_1$. Once the reaction occurs, the molecular $S_3$ will be generated after a fixed time delay $\tau$, and will degrade with the rate $k_2$. 
```math
S_1+S_2 \xrightarrow{k_1}\emptyset,~~\emptyset\stackrel{\tau}\Rightarrow S_3,\\
S_3 \xrightarrow{k_2}\emptyset
```

The species are $S_1,S_2,S_3$. Let $k_1=0.001, k_2 = 0.001，\tau = 0.1.$

```R
tmax <- 150
n_initial <- matrix(c(1000,1000,0),nrow = 3)
t_initial <- 0
S_matrix <- c(-1,-1,0,0,0,-1)
S_matrix <- matrix(S_matrix,nrow = 3) 
S_matrix_delay <- c(0,0,1,0,0,0)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 3)
k <- c(0.001,0.001)
reactant_matrix <- matrix(c(1,1,0,0,0,1),nrow = 3)
delay_type <- matrix(c(2,0),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,0.1)
delaytime_list <- append(delaytime_list,0)
```

Simulate $10^4$ times and calculate the mean value for all the molecule species and probability distribution of $S_3$ at $t = 150$. The number of $S_1$ is the same as the number of $S_2$ at any time.

```R
sample <- 10000
result <- simulation_DelaySSA(algorithm = "DelayMNR", sample_size=sample, tmax=tmax, n_initial=n_initial, t_initial=t_initial, S_matrix=S_matrix, S_matrix_delay=S_matrix_delay, k=k, reactant_matrix=reactant_matrix, delay_type=delay_type , delaytime_list=delaytime_list)
plot_SSA_mean(result = result,t=seq(0, tmax, by = 1) ,n_initial = n_initial,t_initial = 0)
plot_SSA_density(result = result,t_pick = tmax)
```

Check [Tutorials](https://github.com/Zoey-JIN/DelaySSA/blob/main/Tutorials.md) for more details.

## References
[1]Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. The journal of physical chemistry, 81(25), 2340-2361.

[2] Cai, X. (2007). Exact stochastic simulation of coupled chemical reactions with delays. The Journal of chemical physics, 126(12).

[3] Barrio, M., Burrage, K., Leier, A., & Tian, T. (2006). Oscillatory regulation of Hes1: discrete stochastic delay modelling and simulation. PLoS computational biology, 2(9), e117.

[4] Anderson, D. F. (2007). A modified next reaction method for simulating chemical systems with time dependent propensities and delays. The Journal of chemical physics, 127(21).





