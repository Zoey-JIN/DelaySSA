
# Gillespie Algorithm

 Consider a system of $N$ chemical substances with $R$ ongoing chemical reactions, each of which has a corresponding tendency function $f_r(\bm{n})$. The Gillespie algorithm assumes that the time from start to finish for each reaction is negligible. Through random simulations, calculate 1) how much time will pass before the next reaction occurs (i.e. starts and finishes), and 2) which reaction will occur at that future point in time. The following assumptions, sometimes referred to as the basic premises of chemical dynamics, are based on physical principles and serve as the underlying assumptions for methods of simulating chemical reaction systems [1]:
 $$
 f_r(\bm{n}(t)) dt = \text{the probability that  reaction r takes place in a small time interval} ~[t, t + dt)
 $$

Based on this  fundamental assumptions,  $\tau$ and $\mu$ are two independent random variables and following the probability density functions as:
$$
p(\tau|\bm{n},t)=\lambda(\bm{n},t) \exp(-\tau \lambda(\bm{n},t)), ~~\lambda=\sum_{r=1}^{R} f_r(\bm{n},t),
$$
$$
p(\mu|\bm{n},t)=f_r(\bm{n},t)/\lambda(\bm{n},t) .
$$
According this two equations, $\tau$ and $\mu$ can be generated as,
$$
\tau=-\ln(u_1)/\lambda(\bm{n},t) ,
$$
$$
\mu=\text{the integer satisfies $\sum_{r=1}^{\mu-1} f_r(\bm{n},t)< u_2 \lambda(\bm{n},t) \leq \sum_{r=1}^{\mu} f_r(\bm{n},t)$} ,
$$
where $u_1,u_2\sim \text{Uniform}(0,1)$ respectively.

<!-- [Approximate accelerated stochastic simulation of chemically reacting systems]
[Improved leap-size selection for accelerated stochastic simulation] -->



## Algorithm

 1. Initialize. Set $ t = 0 $ and set species number $\bm{n}= n_{initial}$. 

 2. Calculate the propensity function, $f_r$, for each reaction.

 3. Generate two independent, uniform $(0,1)$ random numbers, $u_1$ and $u_2$.

 4. Set $\tau = -\ln(u_1)/\sum_{r=1}^{R} f_r$.

 5. Find the integer $\mu$ which satisfies $\sum_{r=1}^{\mu-1} f_r< u_2 \sum_{r=1}^{R} f_r \leq \sum_{r=1}^{\mu} f_r$/

 6. Set $t = t + \tau$.

 7. Update species number $\bm{n}$ based upon the completion of the reaction $\mu$.

 8. Return to step 2 or quit.




# The Next Reaction Method Algorithm
 
 Let $v_r$, $v'_r\in N^N_{\geq 0}$ be the vectors representing the number of each species consumed and created in the *r*th reaction, respectively. Then, if $N_r(t)$ is the number of initiations of reaction $r$ by time $t$, the state of the system at time $t$ is

$$
\bm{n}(t)=\bm{n}(0)+\sum_{r=1}^R{N_r(t)(v^{'}_r-v_r)}.
$$

 However, based upon the fundamental assumption of stochastic chemical kinetics, $N_r(t)$ is a counting process with intensity $f_r(\bm{n}(t))$ such that $p(N_k(t+\Delta t)-N_k(t)=1)=f_r(\bm{n}(t))\Delta t$, for small $\Delta t$. Therefore, we have

$$
N_r(t)=Y_r\Big(\int^t_0f_r(\bm{n}(s))ds\Big),\tag{1} 
$$

where the $Y_r$ are independent, unit rate Poisson processes. 

Each Poisson process $Y_r$ brings its own time frame. If we define $T_r(t)=\int^t_0f_r(\bm{n}(s))ds$ for each $r$, then it is relevant for us to consider $Y_r(T_r(t))$. We will call $T_r(t)$ the "**internal time**" for reaction $r$.

$\Delta t_r$ notes the gap time which the *r*th reaction needs. $\Delta=\min_r { \Delta t_r }$. For the moment we denote $\overline{t} = t +\Delta $ and the updated propensity functions by $\overline{f_r}$. 

The internal time of the next firing of $Y_r$ has not changed and is still given by $T_r(t) + f_r \Delta t_r$. We also know that the updated internal time of $Y_r$ is given by $T_r(\overline{t}) =T_r(t)+ f_r \Delta $. Therefore, the amount of internal time that must pass before the *r*th reaction fires is given as the difference:
$$
(T_r(t) + f_r \Delta t_r) − (T_r(t)+ f_r \Delta ) = f_r(\Delta t_r − \Delta).
$$
Thus, the amount of absolute time that must pass before the *r*th reaction channel fires,$ \Delta \overline{t}_r$, is given as the solution to $\overline{f_r}\Delta \overline{t}_r = f_r(\Delta t_r − \Delta)$, and we can get:
$$
\overline{\tau_r} = f_r / \overline{f_r}  (\Delta t_r − \Delta) + \overline{t}
    = f_r / \overline{f_r}  ((t+\Delta t_r )− (t+\Delta)) + \overline{t}
    = f_r / \overline{f_r}  (\tau_r − \overline{t}) + \overline{t}
$$

We have therefore found the absolute times of the next firings of reactions $r = µ$ without having to generate any new random numbers.

 Note that after the first timestep is taken in the Next Reaction Method, all subsequent timesteps only demand one random number to be generated.  The Next Reaction Method developed with the notion of a dependency graph and a priority queue in order to increase computational efficiency [2]. This is similar with random numbers needed for each step of the original Gillespie Algorithm.

## Algorithm

 1. Initialize. Set $ t = 0 $ and set species number $\bm{n}= n_{initial}$. 

 2. Calculate the propensity function, $f_r$, for each reaction.

 3. Generate $R$ independent, uniform $(0,1)$ random numbers, $u_r$.

 4. set $\tau_r = -\ln(u_r)/ f_r$.

 5. Set $t = \min_r \{ \tau_r \}$. Here we assume that $\tau_\mu$ is the minimum.

 6. Update species number $\bm{n}$ based upon the completion of the reaction $\mu$.

 7. Recalculate the propensity function, $\overline{f_r}$, for each reaction.

 8. For each $r \neq \mu$, set $\tau_r=(f_r/\overline{f_r})(\tau_r-t)+t$.

 9. For reaction $\mu$, let $u^{'}$ be uniform $(0,1)$ and set $\tau_\mu=-\ln(u^{'})/\overline{f_
 \mu}+t$. If $\tau_r$ is either $NA$ or $Inf$, it also needs to be recalculated in this manner.

 10. For each r, set $f_r = \overline{f_r}$.

 11. Return to step 5 or quit.


# Modified Next Reaction Method Algorithm
According to [3], Modified Next Reaction Method Algorithm that is completely equivalent to Next Reaction Method Algorithm. But  but makes more explicit use of the internal times $T_r$. The main idea of this algorithm is $\Delta t_r = (1/f_r)(P_r − T_r)$.

## Algorithm

 1. Initialize. Set $ t = 0 $ and set species number $\bm{n}= n_{initial}$. For each $r \leq R$, set $P_r = 0$ and $T_r = 0$.

 2. Calculate the propensity function, $f_r$, for each reaction.

 3. Generate $R$ independent, uniform $(0,1)$ random numbers, $u_r$, and set $P_r = -\ln(u_r)$.

 4. Set $ \tau_r = (P_r − T_r)/f_r$.

 5. Set $\tau = \min_r \{  \tau_r \}$. Here we assume that $\tau_\mu$ is the minimum.

 6. Set $t = t + \tau$. And update species number $\bm{n}$ based upon the completion of the reaction $\mu$.

 7. For each r, set $T_r = T_r+f_r\tau$.

 8. For reaction $\mu$, let $u^{'}$ be uniform$(0,1)$ and set $P_\mu = P_\mu - \ln(u^{'})$.

 9. Recalculate the propensity function, $f_r$, for each reaction.

 10. Return to step 5 or quit.

## References

[1] Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. The journal of physical chemistry, 81(25), 2340-2361.

[2] Gibson, M. A., & Bruck, J. (2000). Efficient exact stochastic simulation of chemical systems with many species and many channels. The journal of physical chemistry A, 104(9), 1876-1889.

[3] Anderson, D. F. (2007). A modified next reaction method for simulating chemical systems with time dependent propensities and delays. The Journal of chemical physics, 127(21).


