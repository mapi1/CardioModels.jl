### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ f27d484a-0993-11eb-08c9-030d9f88b754
using Plots, Distributions, PlutoUI, DSP, Statistics, CardioModels

# ╔═╡ 60c79388-12b7-11eb-0493-51151fa697ab
md"""
### About the model

A main feature of the model is that besides parasympathetic action through the baroreflex, sympathetic action is added. Sympathetic activity is triggered if the diastolic pressure falls below a reference value and influences heart rate, inotropy, and systemic resistance with some delay. In this way, a 0.1Hz oscillation is induced.

#### The overall model is formulated as follows

##### Sympathetic activation

$SympD_n = 0 \hspace{5mm} \text{if} \hspace{22mm} D_{n} > D_{\text{ref}}$
$SympD_n = D_{\text{ref}} - D_n \hspace{5mm} \text{if} \hspace{5mm} D_n \leq D_{\text{ref}}$
$Symp_n = \text{atan} \left(G_{SA}(z) \cdot SympD_{n}\right)$


##### Baroreflex

$I_n^{\star} = a_1 S_{n-1} + I_0 + w_I + A_{\rho}^I\sin(2\pi f_{\rho}t + \phi) \hspace{5mm} \text{if} \hspace{5mm}I_{n-1} > 700\ ms$

$I_n^{\star} = a_2 S_{n-2} + I_0 + w_I + A_{\rho}^I\sin(2\pi f_{\rho}t + \phi) \hspace{5mm} \text{if} \hspace{5mm}I_{n-1} \leq 700\ ms$

$I_n = \frac{I_n^{\star}}{1 + \beta Symp_n}$

with $\beta$ the sympathetic 'force' on $I$

##### Windkessel

$D_n = S_{n} \exp\left(\frac{-I_n }{R C (1 + \delta Symp_n)}\right)$

with $RC$ the diastolic runoff and $\delta$ the sympathetic 'force' on $D$

##### Starling/Restitution

$P_n = (\gamma I_n + P_0)(1+ \varepsilon Symp_n) + w_P + A_{\rho}^P\sin(2\pi f_{\rho}t + \phi)$
with $\varepsilon$ the sympathetic 'force' on $P$ (inotropy by symathetics and respiratory BP modulation)
\
The next systole is determined by $S_{n+1} = D_n + P_n$ 
"""

# ╔═╡ 004f4ad1-4160-4e52-8813-70b9c73e3a29
md"# Create a model:"

# ╔═╡ 7efa1b0e-0a5c-11eb-14c9-032c18080c9b
md"""
Diastolic reference value:
\
$(@bind Dref NumberField(30:160, default = 80))
mmHg
"""

# ╔═╡ 25b183da-0998-11eb-0a1f-03534abcdc2e
md"""
Diastolic runoff:
\
$(@bind t_runoff NumberField(1:5000, default = 2000))
ms
"""

# ╔═╡ cdba0c98-120b-11eb-2aa3-d9e84e989773
md"""
Offset:
\
$(@bind offset NumberField(1:100, default = 20))
mmHg
"""

# ╔═╡ ac0c0478-0998-11eb-31bc-e3e4193eb47f
md"""
Baro strength: 
\
a1:
$(@bind a1 NumberField(0:0.01:30, default = 9))
mmHg/ms 
"""

# ╔═╡ f7a5c9ca-0998-11eb-1a7f-c57908a3cee5
md"""
a2:
$(@bind a2 NumberField(0:0.01:30, default = 9))
mmHg/ms
"""

# ╔═╡ 918ea734-120a-11eb-33e7-4340ad35c9f1
md"""
Sympathetic Force on D:
\
$(@bind b_D NumberField(0:0.001:1, default = 0.25))

"""

# ╔═╡ ee91b962-12c3-11eb-24aa-ff59fbe0e53e
md"""
Sympathetic Force on P:
\
$(@bind b_P NumberField(0:0.001:1, default = 0.25))
"""

# ╔═╡ ed8b4472-12c3-11eb-0d55-47162f88dc88
md"""
Sympathetic Force on I:
\
$(@bind b_I NumberField(0:0.001:1, default = 0.25))
"""

# ╔═╡ eb927d16-120a-11eb-2a23-277da77b03f6
md"""
Starling Effect:
\
$(@bind starling NumberField(0:0.001:1, default = 0.02))
mmHg/ms
"""

# ╔═╡ 78a242de-120b-11eb-1662-a948f6a40cf5
w = [0.000,0.00, 0.015, 0.030, 0.02, 0.010];

# ╔═╡ a10eb0f2-50f1-11eb-0de4-f3e26e2ed897
md"""
Breathing period T:
\
$(@bind T_resp NumberField(0:50:10000, default = 4000))
ms
"""

# ╔═╡ 3d310136-0991-11eb-05ea-4351c72e2260
model = KaremakerModel(a1 = a1, a2 = a2, t_runoff = t_runoff, Dref = Dref, w = w, b_D = b_D, b_I = b_I, b_P = b_P, starling = starling, offset = offset, T_resp = T_resp, A_resp = 60, magn_pulse = 5, A_BP = 10);
#model = KaremakerModel();

# ╔═╡ 2bc64aaf-c07a-4e75-84e7-ce71f743f280
md"## Generate data:"

# ╔═╡ 23cc206d-78a3-4be0-b653-a8c5123d8dd9
md"""
N:
$(@bind N NumberField(1:1000, default = 100))
"""

# ╔═╡ 4935f00e-0991-11eb-2017-210b39b3aa56
S, D, P, I, Symp, ρ = predict(model, N, burnIn = 200);

# ╔═╡ b01c33d6-0997-11eb-222c-d34da1f9aa4b
md"""
Schow spectral density:
$(@bind p_psd CheckBox())
"""

# ╔═╡ e3c9adde-0993-11eb-191c-b501deb4a46f
begin
	if !p_psd
		pi = plot(I, title = "RR [ms]", lab = "")
		pbp = plot(S, title = "BP [mmHg]", lab = "", color = :red)
		plot!(pbp, D, lab = "", color = :green)
		plot(pi, pbp, layout = (2,1), dpi = 500)
	else
		psdi = psdplot(I, title = "RR [ms]", method = "welch", lab = "",  fs = 1000/ (mean(I)))
		psdbp = psdplot(S, title = "BP [mmHg]", lab = "", color = :red, fs = 1000/ (mean(I)))
		psdplot!(psdbp, D, lab = "", color = :green, fs = 1000/ (mean(I)))
		plot(psdi, psdbp, layout = (2,1), size = (700, 500), fs = 1000/ (mean(I)))
		vline!([1/(0.001T_resp)], lab = "f_resp")
	end
end

# ╔═╡ 0f3b19aa-120c-11eb-3c2a-6f5e3e1e42fa
md"""
**Experiment on LF / HF ratio as measure for sympatho-vagal balance:**
\
Diastolic runoff as a measure of overall systemic resistence is varied resulting in a change in sympathetic drive, which is activated only over the diastolic pressure in this model. 
"""

# ╔═╡ 4eb77de8-120d-11eb-3d3b-f1414ea32fef
begin
	model2 = KaremakerModel(b_I = 0.5)
	runoffs = range(1250, 2850, length = 100)
	rr,symp_drive, lfhf, sb, db = vagal_balance(55, runoffs, model = model2)
	pv = plot(symp_drive, lfhf, ylab = "LF/HF ratio", lab = "", title = "Sympathovagal Balance")
	pr = plot(symp_drive, rr, ylab = "RR [ms]", lab = "")
	p_BP = plot(symp_drive, sb, lab = "SB")
	plot!(symp_drive, db, lab = "DB", xlab = "sympathetic drive")
	#plot(pv, pr, p_BP, layout = (3,1))
	plot(pv, pr, layout = (2,1))
end

# ╔═╡ Cell order:
# ╠═f27d484a-0993-11eb-08c9-030d9f88b754
# ╟─60c79388-12b7-11eb-0493-51151fa697ab
# ╟─004f4ad1-4160-4e52-8813-70b9c73e3a29
# ╠═3d310136-0991-11eb-05ea-4351c72e2260
# ╟─7efa1b0e-0a5c-11eb-14c9-032c18080c9b
# ╟─25b183da-0998-11eb-0a1f-03534abcdc2e
# ╟─cdba0c98-120b-11eb-2aa3-d9e84e989773
# ╟─ac0c0478-0998-11eb-31bc-e3e4193eb47f
# ╟─f7a5c9ca-0998-11eb-1a7f-c57908a3cee5
# ╟─918ea734-120a-11eb-33e7-4340ad35c9f1
# ╟─ee91b962-12c3-11eb-24aa-ff59fbe0e53e
# ╟─ed8b4472-12c3-11eb-0d55-47162f88dc88
# ╟─eb927d16-120a-11eb-2a23-277da77b03f6
# ╠═78a242de-120b-11eb-1662-a948f6a40cf5
# ╟─a10eb0f2-50f1-11eb-0de4-f3e26e2ed897
# ╟─2bc64aaf-c07a-4e75-84e7-ce71f743f280
# ╟─23cc206d-78a3-4be0-b653-a8c5123d8dd9
# ╠═4935f00e-0991-11eb-2017-210b39b3aa56
# ╟─b01c33d6-0997-11eb-222c-d34da1f9aa4b
# ╟─e3c9adde-0993-11eb-191c-b501deb4a46f
# ╟─0f3b19aa-120c-11eb-3c2a-6f5e3e1e42fa
# ╠═4eb77de8-120d-11eb-3d3b-f1414ea32fef
