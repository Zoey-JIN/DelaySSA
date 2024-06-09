# Basic Concepts
Biological processes encompass a multitude of intricate mechanisms involving various molecules and physical interactions. By conceptualizing these processes as a sequence of discrete chemical reactions, we can mathematically formalize them. This mathematical representation enables precise modeling of the reaction dynamics and facilitates accurate predictions of reactant quantities.

Stochastic Simulation Algorithm (SSA) [1] is a method used to simulate stochastic processes in chemical reaction systems. This algorithm is particularly suitable for systems with a small number of molecules.

However, some chemical reactions, such as gene transcription and translation within living cells, necessitate a defined temporal duration to reach completion following their initiation. As a result, the products of these reactions will manifest only after a time delay which makes the traditional SSA algorithm unsuitable for these reactions. DelaySSA implements a stochastic simulation algorithm (SSA) with delays in R language. It can simulate chemical reaction systems both with and without delays. now DelaySSA supports three exact delay stochastic simulation algorithms, namely, delay direct method `DelayDirect` [2], delay rejection method DelayRejction `DelayReject` [3] and delay modified next reaction method `DelayMNR` [4]. Here we give the basic concepts of these algorithms.

Given a finite set of chemical species $X_i, i = 1, \ldots,N,$ and $R$ chemical reactions, we define the reactions by the notation  

$$
\sum_{i=1}^{N} s_{ir}X_i \xrightarrow{k_r} \sum_{i=1}^{N} s^{'}_{ir}X_i,~~r=1, \ldots,R,
$$  

where $s_{ir}$ and $s'_{ir}$ denote numbers of reactant and product molecules, respectively. $ k_r $ is the reaction rate constant of the $r$-th reacion. And the stoichiometric matrix $S$ is given by
$$
S_{ir}=s^{'}_{ir}-s_{ir},~~r=1, \ldots,R,~~i=1, \ldots,N.
$$

According to [Stochastic processes in physics and chemistry], propensity function $f(\bm{n})$ are in the form of mass-action kinetics type

$$
f_r(\bm{n})_=k_r \Omega \prod_{i=1}^{N} \frac{n_i!}{(n_i-s_{ir})! \Omega^{s_{ir}}},
$$

where $\bm{n} = \left( n_1, \ldots, n_N \right)$, $n_i$ is the number of species $X_i$, $\Omega $ is the volume of a closed compartment.

The time delay could be a fixed number or a stochastic value. According to Barrio et al.[Oscillatory Regulation of Hes1: Discrete Stochastic Delay Modelling and Simulation], reactions with delays are categorized into consuming and nonconsuming reactions. If a delayed reaction is a nonconsuming reactions, it initiates at $t$ and will finish until $t+t_\text{delay}$, then the $\bm{n}$ of the number of species will change only at $t+t_\text{delay}$. If a delayed reaction is a consuming reactions, it initiates at $t$ and will finish until $t+t_\text{delay}$, then the $\bm{n}$ of the number of species will change both at $t$ and $t+t_\text{delay}$.

We can categorize reactions into the following three cases.

Case 1: If reaction $r$ loses the reactant species and gains the product species at the initiation time $t$, we define the reaction r with no delays as ND.

Case 2: If reaction r loses the reactant species and gains the product species at the completion time $t+t_\text{delay}$, we define the reaction r with delays as CD.

Case 3: If reaction r loses the reactant species and gains the product species at the initiation time $t$ and the completion time $t+t_{delay}$, respectively, we define the reaction $r$ with delays as ICD.

## References
[1]Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. The journal of physical chemistry, 81(25), 2340-2361.

[2] Cai, X. (2007). Exact stochastic simulation of coupled chemical reactions with delays. The Journal of chemical physics, 126(12).

[3] Barrio, M., Burrage, K., Leier, A., & Tian, T. (2006). Oscillatory regulation of Hes1: discrete stochastic delay modelling and simulation. PLoS computational biology, 2(9), e117.

[4] Anderson, D. F. (2007). A modified next reaction method for simulating chemical systems with time dependent propensities and delays. The Journal of chemical physics, 127(21).
