# CardioModels.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg) [![codecov.io](http://codecov.io/github/mapi1/CardioModels.jl/coverage.svg?branch=master)](http://codecov.io/github/mapi1/CardioModels.jl?branch=master) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mapi1.github.io/CardioModels.jl/dev)

This little package gives some models of cardiovascular variability, for experimenting and learning.

# Installation

To install this package type the following into the REPL
```julia
using Pkg
Pkg.add("https://github.com/mapi1/CardioModels.jl.git")
```

# Examples

Example notebooks using [Pluto.jl](https://github.com/fonsp/Pluto.jl) are provided in the folder `/examples` showing basic usage examples. The folder contains its own `Project.toml`. By running Julia in the example folder and using

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

all required packages will be installed.


# Models
Below is a minimal introduction of the models, for more details refer to my thesis, the docs or the original publications. The following convention is used for cardiovascular variables:

* I:  RR Intrval series
* S:  Systolic blood pressure
* D:  Diastolic blood pressure
* P:  Pulse Pressure
* <a href="https://www.codecogs.com/eqnedit.php?latex=\rho" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho" title="\rho" /></a>:   Respiration
* T: Time constant as a measure of peripheral resistance R

The relation between different variables is shown in the below figure:
![Notation and relation of cardiovascular variables](https://github.com/mapi1/CardioModels.jl/blob/master/examples/figures/rst.svg)
## DeBoer-model

*Source:* DeBoer, R. W., Karemaker, J. M., & Strackee, J. (1987). Hemodynamic fluctuations and baroreflex sensitivity in humans: a beat-to-beat model. American Journal of Physiology-Heart and Circulatory Physiology, 253(3), H680-H689. ([DOI][deb87])

[deb87]: https://doi.org/10.1152/ajpheart.1987.253.3.H680
## Karemaker-model

*Source:* Karemaker, J. M. (1998). Testing the validity of LF/HF as measure of ‘sympathovagal balance’ in a computer model of cardiovascular control. Proceedings of the IX International Symposium on the Autonomic Nervous System. 
## Baselli-model

*Source:* Baselli, G., Cerutti, S., Civardi, S., Malliani, A., & Pagani, M. (1988). Cardiovascular variability signals: towards the identification of a closed-loop model of the neural control mechanisms. IEEE Transactions on Biomedical Engineering, 35(12), 1033-1046. ([DOI][bas88])

[bas88]: https://doi.org/10.1109/10.8688
