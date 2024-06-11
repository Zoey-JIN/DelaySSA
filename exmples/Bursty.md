# Bursty Model

We study the Bursty model which describes gene expression as occurring in bursts, where multiple mRNA molecules are rapidly synthesized during short periods of high activity and will degrade after the fixed delay time $\tau$. The rate of gene expression occurs is given by the function $f(n)=\frac{ab^n}{(1+b)^{n+1}}$ for any integer $n$. This can be described by the reaction scheme:

$$
\emptyset\stackrel{k_i=\frac{\alpha b^i}{(1+b)^{i+1}}}\longrightarrow iN, ~~iN\stackrel{\tau}\Rightarrow\emptyset,~~i=1,2,3,...
$$

The specie is $N$. Let $i=1,\ldots,30,~a=0.0282,b=3.46,\tau=120$.

```R
j <- 30
tmax <- 200
n_initial <- matrix(c(0),nrow = 1)
t_initial <- 0
S_matrix <- c(1:j)
S_matrix <- matrix(S_matrix,nrow = 1) 
S_matrix_delay <- -c(1:j)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 1)
a <- 0.0282
b <- 3.46
k <- c(sapply(1:j, function(i) a * b^i / (1 + b)^(i + 1)))
product_matrix <- matrix(rep(0,j),nrow = 3)
delay_type <- matrix(rep(c(2),times=j),nrow = 1)
delaytime_list <- list()
for (i in 1:j) {
  delaytime_list <- append(delaytime_list,120) 
}
```

We simulate $10^5$ trajectories and calculate the mean value and probability distribution at $t = 200$.

```R
sample <- 100000
result <- simulation_DelaySSA(algorithm = "DelayMNR", sample_size=sample, tmax=tmax, n_initial=n_initial, t_initial=t_initial, S_matrix=S_matrix, S_matrix_delay=S_matrix_delay, k=k, product_matrix=product_matrix, delay_type=delay_type , delaytime_list=delaytime_list)
```

![Bursty_mean](../figs/Bursty_mean.svg)
![Bursty_density_200s](../figs/Bursty_density_200s.svg)