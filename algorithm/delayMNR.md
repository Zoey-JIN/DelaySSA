# Delay Modified Next Reaction Method Algorithm

 Delay Modified Next Reaction Method Algorithm is modified Next Reaction Method to systems with delays. According to [1], if $T_r$ is the current internal time of $Y_r$, $P_r$ the first internal time after $T_r$ at which $Y_r$ fires, and the propensity function for the *r*th reaction channel is given by $f_r$, then the time until the next initiation of reaction r (assuming no other reactions initiate or complete) is still given by $\Delta t_r = (1/f_r)(P_r − T_r)$.  To each delayed reaction channel we therefore assign a vector, $s_r$, that stores the completion times of that reaction in ascending order. Thus, the time until there is a change in the state of the system, be it an initiation or a completion, will be given by
 $$
 \Delta=\min\{\Delta t_r, s_r[1]-t\}
 $$
 where $t$ is the current time of the system.

## Algorithm

 1. Initialize. Set $ t = 0 $ and set species number $\bm{n}= n_{initial}$. For each $r \leq R$, set $P_r = 0$, $T_r = 0$, and $s_r = [\infty]$.

 2. Calculate the propensity function, $f_r$, for each reaction.

 3. Generate $R$ independent, uniform $(0,1)$ random numbers, $u_r$, and set $P_r = -\ln(u_r)$.

 4. Set $ \tau_r = (P_r − T_r)/f_r$.

 5. Set $\tau = \min_r \{  \tau_r , s_r[1]-t \}$

 6. Set $t = t + \tau$.

 7. If choose the completion of the delayed reaction $\mu$:

      - Update species number $\bm{n}$ based upon the completion of the reaction $\mu$.

      - Delete the first element of $S_\mu$.

 8. Elseif reaction $\mu$ initiated and $\mu\in \text{ND}$

      - Update species number $\bm{n}$ according to reaction $\mu$.

 9. Elseif reaction $\mu$ initiated and $\mu\in \text{CD}$

      - Update $s_\mu$ by inserting $t + \tau_\mu$ into $s_\mu$ ensuring that the times in $s_\mu$ remain in ascending order.

 10. Elseif reaction $\mu$ initiated and $\mu\in \text{ICD}$
    
      - Update species number $\bm{n}$ based upon the initiation of reaction $\mu$.

      - Update $s_\mu$ by inserting $t + \tau_\mu$ into $s_\mu$ ensuring that the times in $s_\mu$ remain in ascending order.

 11. For each r, set $T_r = T_r + f_r \tau$.

 12. If reaction $\mu$ initiated, let $u^{'}$ be uniform$(0,1)$ and set $P_\mu = P_\mu - \ln(u^{'})$.

 13. Recalculate the propensity functions, $f_r$.

 14. Return to step 4 or quit.

## References
[1] Anderson, D. F. (2007). A modified next reaction method for simulating chemical systems with time dependent propensities and delays. The Journal of chemical physics, 127(21).
