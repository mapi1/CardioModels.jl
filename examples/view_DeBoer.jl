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

# ╔═╡ 5c3096b0-01b6-11eb-2c01-5ba0aac27b70
using CardioModels, PlutoUI, Plots, Distributions

# ╔═╡ 0c720f5f-a8be-4d8e-a803-461118c3b2be
md"
# DeBoer Model
"

# ╔═╡ 2e8472f6-01b4-11eb-1f49-b35e9ee7c70b
model = DeBoerModel();

# ╔═╡ 373eadc4-01b6-11eb-0ef3-73d0e4518a24
md"""
N:
$(@bind N NumberField(1:1000, default = 100))
"""

# ╔═╡ 232d91d8-01b6-11eb-08d6-3343be95b77f
S, D, I, T, ρ = predict!(model, N, burnIn = 100);

# ╔═╡ b19bf544-b153-41a0-a013-992da6cdd0cd
md"""
Show Spectrum:
$(@bind psd CheckBox())
"""

# ╔═╡ 8eb6f0f6-0a4b-4fff-bf25-02a37be7b892
md"""
DPI:
$(@bind dpi NumberField(1:1000, default = 500)) adjust if plots are to small/ big.
"""

# ╔═╡ 579f32aa-01b6-11eb-3a79-8f2207aeb03c
begin
	if psd
		println("not here yet!")
	else
		s = plot(S, lab = "", ylabel = "mmHg", color = :red)
		plot!(D, lab = "", title = "BP", ylabel = "mmHg", color = :green)
		t = plot(T, lab = "", title = "T", ylabel = "a.u.", color = :orange)
		i = plot(I, lab = "", title = "I", ylabel = "ms", color = :blue)
		r = plot(ρ, lab = "", title = "ρ", ylabel = "a.u.", xlab = "beats", color = :black)
		plot(s,t,i,r, layout = (4,1), dpi = dpi) # Change dpi if the plot is to big/small
	end
end

# ╔═╡ ac3a7eac-6bba-11eb-37ec-49854b22053f
md"
# Simulating a phenylephrine injection

"

# ╔═╡ 53b158a2-6bbb-11eb-0fa7-87390c33cb64
begin
	model_phe = DeBoerModel(a0 = 9, A = 1,  noiseI = Normal(0, 25), noiseS = Normal(0, 2))
	# pre-injection
	Npre = 50
	S1,D1,I1,T1,R1 = predict!(model_phe, Npre, burnIn = 100)
	
	# injection
	duration = 15
	S2,D2,I2,T2,R2 = phenylephrine(model_phe, 1000, duration)
	
	# post-injection
	Npost = 10 # length of plateau
	S3,D3,I3,T3,R3 = predict!(model_phe, Npost, burnIn = 0)
	
	# back to default
	S4,D4,I4,T4,R4 = phenylephrine(model_phe, -1000, 10)
	
	plot([S1; S2; S3; S4], lab = "", title = "BP", ylab = "mmHg", color = :red)
	pBP = plot!([D1; D2; D3; D4], lab = "", color = :green)
	pIBI = plot([I1; I2; I3; I4], lab = "", ylab = "ms", title = "I", color = :blue)
	plot(pBP, pIBI, layout = (2,1))
	vline!([Npre Npre], color = :black, lab = "", dpi = dpi)
end

# ╔═╡ b312bc92-6bbf-11eb-131e-c946cf3eb23f
begin
	xbp = 0:30
	yms = 9 .* xbp
	scatter(S2 .- minimum(S2), I2 .- minimum(I2), lab = "", title = "Simulated response")
	plot!(xbp, yms, color = :black, lab = "BRS 9 ms/mmHg", ylab = "ΔI [ms]", xlab = "ΔS [mmHg]", dpi = dpi)
end

# ╔═╡ Cell order:
# ╠═5c3096b0-01b6-11eb-2c01-5ba0aac27b70
# ╟─0c720f5f-a8be-4d8e-a803-461118c3b2be
# ╠═2e8472f6-01b4-11eb-1f49-b35e9ee7c70b
# ╟─373eadc4-01b6-11eb-0ef3-73d0e4518a24
# ╠═232d91d8-01b6-11eb-08d6-3343be95b77f
# ╟─b19bf544-b153-41a0-a013-992da6cdd0cd
# ╟─8eb6f0f6-0a4b-4fff-bf25-02a37be7b892
# ╠═579f32aa-01b6-11eb-3a79-8f2207aeb03c
# ╠═ac3a7eac-6bba-11eb-37ec-49854b22053f
# ╠═53b158a2-6bbb-11eb-0fa7-87390c33cb64
# ╠═b312bc92-6bbf-11eb-131e-c946cf3eb23f
