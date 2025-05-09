# InjectivityStiefel.jl
This repository is associated to the publication ["P.-A. Absil, S. Mataigne, The Ultimate Upper Bound on the Injectivity Radius of the Stiefel Manifold, SIAM Journal on Matrix Analysis and Applications 46(2), 1145-1167, 2025"](https://doi.org/10.1137/24M1644808).
## Installation
In order to use this repository, the user should have a Julia installation where the following packages are installed: `LinearAlgebra`, `SkewLinearAlgebra`, `Plots`, `LaTeXStrings`, `Colors`. These packages are easily obtained from the package installation environment as follows. In Julia REPL, press `]` to access the installation environment and for each package, do
```julia
(@v1.6) pkg> add Name_of_Package
```
## Use
This repository contains de following folders: `src`, `plots` and `figures`.  
* The folder `src` contains the algorithms described in the paper, namely Algorithm 8.1 under the name `checkradius`. These are the routines needed to produce the results.
  * `skewlog.jl` contains a routine to compute the skew-symmetric matrix logarithm of an orthogonal matrix.
  * `check_injectivity_radius.jl` contains the algorithm corroborating the injectivity radius on the Stiefel manifold $\mathrm{St}(n,p)$ and subroutines (Algorithm 8.1).
* The folder `plots` contains files to perform numerical experiments with the algorithms from the folder `src`.
  * `plot_inj_beta.jl` runs an experiment on $\mathrm{St}(n,p)$ using the routine `check_radius`. The routine is called for various values of parameter $\beta>0$, where $\beta$ is the parameter of the metric on   $\mathrm{St}(n,p)$ (see paper). For each $\beta$, the value $\rho$ is set to the conjectured injectivity radius $î_\beta$ (see Theorem 7.1) and $î_\beta+\delta$. For any $\delta>0$, `check_radius` should stop in a finite number of iteration with probability $1$ (see Section 8). When $\rho=î_\beta$, `check_radius` should never return if Conjecture 8.1 is true (i.e., should reach `itermax`). `plot_inj_beta.jl`  returns a plot with, for each pair $(\beta,\rho)$, a white dot if the algorithm returned in finite time and a black dot if the algorithm reached `itermax`.
  * `plot_iteration_count.jl` also runs an experiment on $\mathrm{St}(n,p)$ using the routine `check_radius`. Given $\beta>0$, and various values $\delta>0$, `plot_iteration_count.jl` displays how many iterations are needed to return when `check_radius` is called with $\rho=î_\beta+\delta$. As $\delta\rightarrow 0$, we observe that the number of iterations needed tends to infinity.
* The folder `figures` contains output figures of the file `plot_inj_beta.jl` set with $(n,p) \in \{(4,2),(4,3),(5,3)\}$. On all figures, `check_radius` reached itermax for $\rho = î_\beta$ --- that is, could never contradict the conjecture $î_\beta=\mathrm{inj}_{\mathrm{St}(n,p)}$.

## Bibtex
If you use the content of this repository, please cite:
```
@article{doi:10.1137/24M1644808,
author = {Absil, P.-A. and Mataigne, Simon},
title = {The Ultimate Upper Bound on the Injectivity Radius of the Stiefel Manifold},
journal = {SIAM Journal on Matrix Analysis and Applications},
volume = {46},
number = {2},
pages = {1145-1167},
year = {2025},
doi = {10.1137/24M1644808},
URL = {https://doi.org/10.1137/24M1644808},
eprint = {https://doi.org/10.1137/24M1644808}
}
```

