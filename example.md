# Oscillation model
In this section, we will give an example to illustrate the use of DelaySSA for defining chemical reaction models, solving associated problems, and visualizing the results. The source code for all examples and figures can be found in the GitHub examples folder. We focus on a specific case known as the Oscillation model, which involves the following steps: (i) a protein $X$ is transcribed by a promoter, (ii) after a fixed time delay $\tau$, $X$ is converted into protein $Y$ through a series of unspecified biochemical processes, and (iii) $Y$ then binds to the promoter, reducing the transcription rate of $X$. This process can be represented by the following reaction scheme:
$$
\begin{aligned}
\emptyset \xrightarrow{J_1(Y)} X,
X\stackrel{\tau}\Rightarrow Y,
Y\xrightarrow{J_2(Y)} \emptyset,
\end{aligned}
$$
The function $J_1(Y)$ and $J_2(Y)$ is defined as follows:
$$
\begin{aligned}
J_1(Y)=k_1S\frac{K^p_d}{K^p_d+Y^p},\\
J_2(Y)=k_2E_T\frac{Y}{K_m+Y}.
\end{aligned}
$$

# Model


## Non-Delay part
```
product_list <- matrix(c(0,0,0,1),nrow = 2)
fun_fr <- function(k,n){
  k <- c(1/(1 + (n[2])^2), 1/(1 + n[2]))
  return(k*propensity(n,product_list))
}
S_matrix <- c(1,0,0,-1)
S_matrix <- matrix(S_matrix,nrow = 2)
```

Then we define initial conditions
```
tmax <- 400
n_initial <- matrix(c(0,0),nrow = 2)
t_initial <- 0
```
## Delay part
```
S_matrix_delay <- c(-1,1,0,0)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 2)
delay_relative <- list(c(2,0),c(20,0))
```

define the reaction
```
sample = 1000
reaction <- function(sample, tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, fun_fr, delay_relative){
  result <- lapply(1:sample, function(x) simulate_reaction_delay_direct(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, fun_fr, delay_relative))
  return(result)
}
```

Finally simulate the reaction and plot the results

# Bursty model

We study the Bursty model which describes gene expression as occurring in bursts, where multiple mRNA molecules are rapidly synthesized during short periods of high activity and will degrade after the fixed delay time $\tau$. The rate of gene expression occurs is given by the function $f(n)=\frac{ab^n}{(1+b)^{n+1}}$ for any integer $n$. This can be described by the reaction scheme:

$$
\begin{aligned}
&\emptyset\stackrel{\frac{\alpha b^i}{(1+b)^{i+1}}}\longrightarrow iN,i=1,2,3,...\\ &N\stackrel{\tau}\Rightarrow\emptyset.
\end{aligned}
$$

# Refractory model
We also study the refractory model, which was devised to explain the experimental observation that the distribution of ‘‘off’’ intervals is not exponential but rather has a peak at a nonzero value. The gene state can change between $G_0$, $G_1$ and $G_2$, but gene expression only occurs at the state of $G_2$. The the mRNA will degrade after the fixed delay time $\tau$. This reaction scheme can be illustrated by
$$
\begin{aligned}
G_0\xrightarrow{\sigma_0} G1,~
G_1\xrightarrow{\sigma_1}G_2,~
G_2\xrightarrow{\sigma_2}G_0,~
G_2\xrightarrow{\rho}G_2+N,~
N\stackrel{\tau}\Rightarrow\emptyset.
\end{aligned}
$$

# Birth Death model
We consider a non-Markovian system where mRNA are transcript at a rate $\rho$ and are removed from the system (degraded) after a fixed time delay $\tau$. For each mRNA, its $\tau$ follow a gamma distribution. This reaction scheme can be illustrated by
$$
\emptyset\stackrel{\rho}\rightarrow N, 
N\stackrel{\tau}\Rightarrow\emptyset.
$$

where $\tau\sim\text{Gamma}(7,1)$

# A reaction with two channels
感觉这个可以当做tutorial里面解释代码的例子

We study the reaction with two non-delay channels and one delay channel. This model describes that molecular $S_1$ binds $S_2$ and then disappear with the reaction rate $c_1$. Once the reaction occurs, the molecular $S_3$ will be generated after a fixed time delay $\tau$, and will degrade with the rate $c_2$. This procedure can be described by

$$
S_1+S_2 \xrightarrow{c_1}\emptyset\\
\emptyset\stackrel{\tau}\Rightarrow S_3\\
S_3 \xrightarrow{c_2}\emptyset
$$







