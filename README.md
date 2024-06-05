# DelaySSA

DelaySSA implements a stochastic simulation algorithm (SSA) with delays in R. It can simulate chemical reaction systems both with and without delays.


## Install
Download the installation package from GitHub
```
devtools::install_github("Zoey-JIN/DelaySSA")
```
Or download and install locally
```
devtools::install("~/DelaySSA")
```
Then load
```
library("DelaySSA")
```
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
