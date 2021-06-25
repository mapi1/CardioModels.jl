# DeBoer Model
DeBoer, R. W., Karemaker, J. M., & Strackee, J. (1987). Hemodynamic fluctuations and baroreflex sensitivity in humans: a beat-to-beat model. American Journal of Physiology-Heart and Circulatory Physiology, 253(3), H680-H689. ([DOI](https://doi.org/10.1152/ajpheart.1987.253.3.H680))

## Theory

#### Effective pressure:
The systolic pressure is transformed to be in accordance with the sigmoid-shaped activation curve of the baroreflex. The constant $S_0$ defines the working point.

$S_{n}^{\star} = f(S_n)\hspace{5mm}\text{with} \hspace{5mm} f(S_n) = S_0 + 18\text{atan}\left(\frac{S_{n} -S_0 }{18}\right)$

#### Baroreflex on Heart Rate
The transfer function $G_a(z)$ relates the effective systolic pressure to the RR interval. It models fast vagally mediated and slower, sympathetic, baroreflex effects.

$I_{n} = G_a (z) S_{n}^{\star}$

#### Baroreflex on Peripheral Resistance
The peripheral resistance is influenced by the baroreflex mediated through the sympathetic nervous system reflected in the transfer function $G_b(z)$. $T_0$ is an empirical offset constant (around 3000 ms).

$T_{ n } = R_n C = T_{ 0 } + G_b(z)S_{n}^{\star}$

#### Properties of Myocardium & Respiration
Starling's law governs the pulse pressure through the coupling constant $\gamma$. Respiration also enters through the pulse pressure, modeled as a sine. 

$P_{ n } = \gamma I_{n}  + A_{\rho}\sin(2\pi f_{\rho} t) \hspace{5mm} \text{with} \hspace{5mm} t = \sum I_n$

#### Windkessel
Diastolic pressure is derived from the two-element Windkessel model.

$D_{ n } = S_{n} \exp \left( - \frac{I_{n}}{ T_{n} } \right)$

#### Time Propagation 
The next beat is calculated through the standard relation:

$S_{n+1} = P_{n} + D_{n}$


## Usage

```@docs
DeBoerModel
```

```@docs
predict(::DeBoerModel, ::Int)
```

```@docs
predict!(::DeBoerModel, ::Int)
```