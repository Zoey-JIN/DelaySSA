# Basic Concepts
Biological processes encompass a multitude of intricate mechanisms involving various molecules and physical interactions. By conceptualizing these processes as a sequence of discrete chemical reactions, we can mathematically formalize them. This mathematical representation enables precise modeling of the reaction dynamics and facilitates accurate predictions of reactant quantities.

Stochastic Simulation Algorithm (SSA) [[1]Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. The journal of physical chemistry, 81(25), 2340-2361.] is a method used to simulate stochastic processes in chemical reaction systems. This algorithm is particularly suitable for systems with a small number of molecules.

However, some chemical reactions, such as gene transcription and translation in living cells, require a certain amount of time to complete after initiation. Thus, the products of these reactions will appear after a delay.

DelaySSA implements a stochastic simulation algorithm (SSA) with delays in R. It can simulate chemical reaction systems both with and without delays. Gillespie’s exact stochastic simulation algorithm has been widely used to simulate the stochastic dynamics of chemically reacting systems.

Based on the same fundamental premise of stochastic kinetics used by Gillespie, an exact SSA for chemical reaction systems with delays called the direct method is developed by Cai[Exact stochastic simulation of coupled chemical reactions with delays]. In his paper, it points out that an algorithm modified from Gillespie’s SSA by Barrio et al.[Oscillatory Regulation of Hes1: Discrete Stochastic Delay Modelling and Simulation] is also an exact SSA  called the rejection method for chemical reaction systems with delays. But the rejection method requires more random variables than the direct method. And another algorithm is extended from the modified Next Reaction Method to systems known as the Next Reaction Method for systems with delays.[A modified Next Reaction Method for simulating chemical systems with time dependent propensities and delays]

Given a finite set of chemical species $X_i, i = 1, \ldots,N,$ and R chemical reactions, we define the reactions by the notation  

$$
\sum_{i=1}^{N} s_{ir}X_i \xrightarrow{k_r} \sum_{i=1}^{N} s^{'}_{ir}X_i,~~r=1, \ldots,R.
$$  

$s_{ir}$ and $s'_{ir}$ denote numbers of reactant and product molecules, respectively. $ k_r $ is the reaction rate constant of the $r$-th reacion. And the stoichiometric matrix $S$ is given by
$$
S_{ir}=s^{'}_{ir}-s_{ir},~~r=1, \ldots,R,~~i=1, \ldots,N.
$$

Propensity functions are in the form of mass-action kinetics type
[Stochastic processes in physics and chemistry]
$$
f_r(\bm{n})_=k_r \Omega \prod_{i=1}^{N} \frac{n_i!}{(n_i-s_{ir})! \Omega^{s_{ir}}}
$$

where $\bm{n} = \left( n_1, \ldots, n_N \right)$, $n_i$ is the number of species $X_i$, $\Omega $ is the volume of a closed compartment.

Some reactions are expected to include delays. Such as gene transcription and translation, these chemical reactions may require a specific duration to complete after initiation. And the products of such reactions will emerge after distinct delays.

The delay could be a fixed number or a stochastic value. According to Barrio et al.[Oscillatory Regulation of Hes1: Discrete Stochastic Delay Modelling and Simulation], reactions with delays are categorized into consuming and nonconsuming reactions. If a delayed reaction is a nonconsuming reactions, it initiates at $t$ and will finish until $t+t_\text{delay}$, then the $\bm{n}$ of the number of species will change only at $t+t_\text{delay}$. If a delayed reaction is a consuming reactions, it initiates at $t$ and will finish until $t+t_\text{delay}$, then the $\bm{n}$ of the number of species will change both at $t$ and $t+t_\text{delay}$.

We can categorize reactions into the following three cases.

Case 1: If reaction r loses the reactant species and gains the product species at the initiation time $t$, we define the reaction r with no delays as ND.

Case 2: If reaction r loses the reactant species and gains the product species at the completion time $t+t_\text{delay}$, we define the reaction r with delays as CD.

Case 3: If reaction r loses the reactant species and gains the product species at the initiation time $t$ and the completion time $t+t_{delay}$, respectively, we define the reaction r with delays as ICD.

## References
[1] Cai, X. (2007). Exact stochastic simulation of coupled chemical reactions with delays. The Journal of chemical physics, 126(12).

[2] Barrio, M., Burrage, K., Leier, A., & Tian, T. (2006). Oscillatory regulation of Hes1: discrete stochastic delay modelling and simulation. PLoS computational biology, 2(9), e117.



