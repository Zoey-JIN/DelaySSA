# Delay Rejection Method Algorithm

 Based upon the no delay algorithm, we see that simulation methods for systems with delays need to calculate when reactions initiate and store when they complete. However, because of the delayed reactions, the propensity functions can change between initiation times. Bratsun et al.[Delay-induced stochastic oscillations in gene regulation] and Barrio et al[Oscillatory regulation of Hes1: discrete stochastic delay modelling and simulation] used an algorithm for computing the initiation times that is exactly like the original Gillespie Algorithm except that if there is a stored delayed reaction set to finish within a computed timestep, then the computed timestep is discarded, and the system is updated to incorporate the stored delayed reaction. The algorithm then attempts another step starting at its new state. This algorithm will be referred to as the Rejection Method.


## Algorithm

 1. Initialize. Set $ t = 0 $ and set species number $\bm{n}= n_{initial}$.

 2. Calculate propensity functions $f_r(t), r=1, \ldots,R$. 

 3. Generate $u_1$ from a uniform(0,1) random variable, and set $\tau = -\ln(u_1)/\sum_{r=1}^{R} f_r$.

 4. If there is a delayed reaction to finish in $[t,t+\tau)$
    - Discard $\tau$.
    - Update $t$ to be the time of the first next delayed reaction and update the species number.
    - Return to step 2 or quit

 5. Else if there is no delayed reaction in $[t,t+\tau)$. 
    - Generate $u_2$ from a uniform(0,1) random variable, and find the integer $\mu$ which satisfies $\sum_{r=1}^{\mu-1} f_r< u_2 \sum_{r=1}^{R} f_r \leq \sum_{r=1}^{\mu} f_r$.

 6. Update according to the type of reaction $\mu$ belongs to: if the reaction $\mu$ belongs to type ND, update species number $\bm{n}$; if the reaction belongs to type CD, store the time $t+\tau_\mu$; if the reaction belongs to type ICD, update species number $\bm{n}$ and store the time $t+\tau_\mu$.

 7. Set the time $t= t+\tau$.

 8. Return to step 2 or quit.