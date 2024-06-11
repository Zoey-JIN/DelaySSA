# Delay Direct Method Algorithm

Consider that $N_d$ delay reactions are ongoing at the time $t$. The delay reactions will complete at $t+T_1,\ldots,t+T_{N_d}$, where $T_1 \leq T_2,\ldots,t+T_{N_d}$. As in the derivation of Gillespie’s exact Stochastic Simulation Algorithm (SSA), $p(\tau, \nu) d\tau $ can be found from the fundamental assumption as $p(\tau, \nu) d\tau = p_0(\tau) f_\mu(t + \tau) d\tau$, where $p_0(\tau)$ is the probability that no reaction will happen in the time interval $[t,t+\tau)$. The delay effects the propensity function. So $p_0(\tau)$ comes to $ \exp(-\sum_{j=0}^{i-1}\lambda(t+T_j)(T_{j+1}-T_j)-\lambda(t+T_i)(\tau-T_{i})),~~\tau \in [T_i, T_{i+1}), ~~i=0,\ldots,N_d$, where the exponent assume equal to zero when $i=0$.
$$
p(\tau|\bm{n},t)=\lambda(t + T_i) \exp ( - \sum_{j=0}^{i-1} \lambda(t + T_j)(T_{j+1} - T_j) - \lambda(t + T_i)(\tau - T_i) ),~~\lambda(t + T_i)=\sum_{r=1}^{R} f_r(t + T_i),
$$
$$
p(\mu|\bm{n},t)=f_r(t + T_i)/\lambda(t + T_i),~~\mu= 1, \ldots,R, ~~\tau \in [T_i, T_{i+1}), ~~i=0,\ldots,N_d.
$$
According this two equations, $\tau$ and $\mu$ can be generated as,
$$
\tau=T_i+\frac{-\ln(1-u_1)-\sum_{j=0}^{i-1} \lambda(t + T_j)(T_{j+1} - T_j)}{\lambda(t + T_i)} ,
$$
$$
\mu=\text{the integer satisfies $\sum_{r=1}^{\mu-1} f_r(t + T_i)< u_2 \lambda (t + T_i) \leq \sum_{r=1}^{\mu} f_r(t + T_i)$} ,
$$
where $u_1,u_2\sim \text{Uniform}(0,1)$ respectively.

## Algorithm

Suppose that at time $t$ there are ongoing delayed reactions set to complete at times $t+T_1, t+T_2, \ldots, t+T_d$. Define $T_0=0$ and $T_{d+1}=\infty$.

Define *Tstruct*, whose $i$-th $(i=1,\dots,d)$ row stores $T_i$ and the index, $\mu_i$, of the reaction that $T_i$ is associated with.

 1. Initialize. Set $ t = 0 $ and set species number $n = n_\text{initial}$. Create a empty *Tstruct*.

 2. Calculate propensity functions $f_r(t), r=1, \ldots,R$. 

 3. Generate  $\tau$.
    
      - Generate an independent $\text{Uniform}(0,1)$ random number $u_1$

      - If *Tstruct* is empty, it means there is no ongoing delayed reaction
          - set $\tau = -\ln(u_1)/\sum_{r=1}^{R} f_r$.

      - Else

          - set $a_{mask}=0$

          - Set $i=0$, $F=0$ and $a_t = \sum_{r=1}^{R} f_rT_1$.
        
          - While $F < u_1$
              
              - Calculate $F=1-e^{-a_t},i=i+1$.

              - Calculate propensity $f_r(t+T_i)$ due to the finish of the delayed reaction at $t+T_{i}$ and calculate $\sum_{r=1}^{R}f_r(t+T_i)$.           

              - Set $a_{mask}=a_t$. Update $a_t=a_t+\sum_{r=1}^{R}f_r(t+T_i)(T_{i+1}-T_i)$.

              - If $i>1$, update species number $n$ due to the delay reaction at $t+T_{i-1}$

          - EndWhile

          - Set $i=i-1$.

          - Calculate $\tau=T_i-(\ln(1-u_1)+a_\text{mask}-\sum_{r=1}^{R}f_r(t+T_i)(T_{i+1}-T_i))/\sum_{r=1}^{R}f_r(t+T_i)$.

      - EndIf

 4. If $\tau\in[T_i,T_{i+1})$, delete the columns $1,\ldots,i$ of $T_i$ and set $T_j=T_j-\tau$.

 5. Generate $u_2$ from a uniform(0,1) random variable, and find the integer $\mu$ which satisfies $\sum_{r=1}^{\mu-1} f_r< u_2 \sum_{r=1}^{R} f_r \leq \sum_{r=1}^{\mu} f_r$.

 6. Update according to the type of reaction $\mu$ belongs to: if the reaction $\mu$ belongs to type ND, update species number $\bm{n}$; if the reaction belongs to type CD, store the time $t+\tau_\mu$; if the reaction belongs to type ICD, update species number $\bm{n}$ and store the time $t+\tau_\mu$. If it is a delay reaction, insert $\tau_\mu$ and the reaction $\mu$ into *Tstruct*, ensuring that the times in *Tstruct* remain in ascending order.

 7. Set $t=t+\tau$.

 8. Return to Step 3 or quit.

Remark. Notice that in the above pseudocode, we modified the Step 4 in the original algorithm for computational efficiency, but both are equivalent. More details are illustrated in [1].

## References
[1] Cai, X. (2007). Exact stochastic simulation of coupled chemical reactions with delays. The Journal of chemical physics, 126(12).




