# Baselli Model
Baselli, G., Cerutti, S., Civardi, S., Malliani, A., & Pagani, M. (1988). Cardiovascular variability signals: towards the identification of a closed-loop model of the neural control mechanisms. IEEE Transactions on Biomedical Engineering, 35(12), 1033-1046. ([DOI](https://doi.org/10.1109/10.8688)))
 

## Theory
The model is defined through the following three equations, where $G$ denotes an all-zero and $M$ an all-pole transfer function:


```math
\begin{aligned}
S_n &= G_{SS}(z)S_n \quad + \quad G_{SI}(z)I_n \quad  + \quad G_{S\rho}(z) \rho_n  \quad  + \quad  M_{S}(z) w_{S} \\
I_n &= G_{IS}(z)S_n \quad  + \hspace{29mm}  G_{I\rho}(z) \rho_n \quad  + \quad  M_{I}(z) w_{I}\\
\rho_n &= \hspace{76mm}M_{\rho}(z) w_{\rho}
\end{aligned}
```

It can be identified from data as an ARXAR-model using the generalized least squares method (via [arxar](https://github.com/baggepinnen/ControlSystemIdentification.jl)).

## Usage

```@docs
BaselliModel
```
The model can be plotted:
```@setup 1
using Plots, CardioModels
```

```@example 1
model = getModel(1)
plot(model)
```

Predefined models can be obtained with:
```@docs
getModel
```

```@docs
predict(::BaselliModel, ::Int)
```

```@docs
fitBaselli
```

```@docs
postprocess
```

