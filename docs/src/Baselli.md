# Baselli Model
The model is defined through the following three equations, where $G$ denotes an all-zero and $M$ an all-pole transfer function:

$S_n = G_{S S}(z) S_n \quad + \quad G_{SI}(z)I_n \quad  + \quad G_{S\rho}(z) \rho_n  \quad  + \quad  M_{S}(z) w_{S}$
$I_n = G_{IS}(z)S_n \quad  + \hspace{40mm}  G_{I\rho}(z) \rho_n \quad  + \quad  M_{I}(z) w_{I}$
$\rho_n = \hspace{104mm}M_{\rho}(z) w_{\rho}$

It can be identified from data as an ARXAR-model using the generalized least squares method.