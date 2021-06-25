# Karemaker Model

Karemaker, J. M. (1998). Testing the validity of LF/HF as measure of ‘sympathovagal balance’ in a computer model of cardiovascular control. Proceedings of the IX International Symposium on the Autonomic Nervous System.

## Theory
A main feature of the model is that besides parasympathetic action through the baroreflex, sympathetic action is added. Sympathetic activity is triggered if the diastolic pressure falls below a reference value and influences heart rate, inotropy, and systemic resistance with some delay. In this way, a 0.1Hz oscillation is induced.

#### Sympathetic activation

```math
\begin{aligned}
SympD_n &= 0 \hspace{18mm} \text{if} \hspace{5mm} D_{n} > D_{\text{ref}} \\
SympD_n &= D_{\text{ref}} - D_n \hspace{5mm} \text{if} \hspace{5mm} D_n \leq D_{\text{ref}}\\
Symp_n &= \text{atan} \left(G_{SA}(z) \cdot SympD_{n}\right)
\end{aligned}
```



#### Baroreflex

$I_n^{\star} = a_1 S_{n-1} + I_0 + w_I + A_{\rho}^I\sin(2\pi f_{\rho}t + \phi) \hspace{5mm} \text{if} \hspace{5mm}I_{n-1} > 700\ ms$

$I_n^{\star} = a_2 S_{n-2} + I_0 + w_I + A_{\rho}^I\sin(2\pi f_{\rho}t + \phi) \hspace{5mm} \text{if} \hspace{5mm}I_{n-1} \leq 700\ ms$

$I_n = \frac{I_n^{\star}}{1 + \beta Symp_n}$

with $\beta$ the sympathetic 'force' on $I$

#####
$D_n = S_{n} \exp\left(\frac{-I_n }{R C (1 + \delta Symp_n)}\right)$

with $RC$ the diastolic runoff and $\delta$ the sympathetic 'force' on $D$

#### Starling/Restitution

$P_n = (\gamma I_n + P_0)(1+ \varepsilon Symp_n) + w_P + A_{\rho}^P\sin(2\pi f_{\rho}t + \phi)$
with $\varepsilon$ the sympathetic 'force' on $P$ (inotropy by symathetics and respiratory BP modulation)
\
The next systole is determined by $S_{n+1} = D_n + P_n$ 

## Usage

```@docs
KaremakerModel
```

```@docs
predict(::KaremakerModel, ::Int)
```