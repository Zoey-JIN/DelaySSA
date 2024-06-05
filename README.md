# DelaySSA

brief illustrasion of this R package DelaySSA

DelaySSA implements a stochastic simulation algorithm (SSA) with delays in R. It can simulate chemical reaction systems both with and without delays.
DelaySSA用r语言的方式实现了对与一些有延迟的反应的ssa模拟。针对这些涉及或者不涉及的反应用包的方式进行模拟。


## Install

how to install 

## Examples

### Telegraph Model

A brief explanation of telegraph model

$$
G\xrightleftharpoons[\sigma_\text{on}]{\sigma_\text{off}}G^\star,~~G\xrightarrow{\rho}G+N,~~N\xRightarrow[\tau]{}\varnothing\tag{1}
$$

distributions plot

### Oscillation Model

A brief explanation of telegraph model

$$
\varnothing \xrightarrow{J_1(Y)} X,~~X \stackrel{\tau}{\Rightarrow} Y,~~ Y\xrightarrow{J_2(Y)} \varnothing
$$

where $J_1(Y)=k_1S K^p_d/(K^p_d+Y^p)$, $J_2(Y)=k_2E_T Y/(K_m+Y)$.

mean value plot

Check [Tutorials](https://github.com/Zoey-JIN/DelaySSA/blob/main/Tutorials.md) for more details.
## References
[1]

[2]...
