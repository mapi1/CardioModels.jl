### A Pluto.jl notebook ###
# v0.14.1

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
using PlutoUI, Plots, KardioUtils, MonteCarloMeasurements, ControlSystems

# ╔═╡ 53b158a2-6bbb-11eb-0fa7-87390c33cb64
begin
	using Distributions
	model_phe = DeBoerModel(a0 = 9, A = 1,  noiseI = Normal(0, 25), noiseS = Normal(0, 2))
	# pre
	Npre = 100
	S1,D1,I1,T1,R1 = predict!(model_phe, Npre, burnIn = 500)
	
	# injection
	S2,D2,I2,T2,R2 = phenylephrine(model_phe, 1000, 10)
	
	# post
	Npost = 20
	S3,D3,I3,T3,R3 = predict!(model_phe, Npost, burnIn = 200)
	
	plot([S1; S2; S3], lab = "", title = "BP", ylab = "mmHg", color = :red)
	pBP = plot!([D1; D2; D3], lab = "", color = :red)
	pIBI = plot([I1; I2; I3], lab = "", ylab = "ms", title = "RR Interval")
	plot(pBP, pIBI, layout = (2,1))
	vline!([Npre Npre], color = :black, lab = "")
end

# ╔═╡ a96c0a2a-01b3-11eb-3975-5953fd76dd7b
include("./Karemaker-DeBoer/DeBoer.jl");

# ╔═╡ 02590918-6c43-11eb-3f54-ab5c83dd1171
include("./Tools/pdc.jl")

# ╔═╡ 16771e7c-6c42-11eb-2b63-6f23d0ffa8d7
cd("../")

# ╔═╡ 0c720f5f-a8be-4d8e-a803-461118c3b2be
md"
# DeBoer Model
Modelling cardiorespiratory varibilities and propagating uncertainty including:
"

# ╔═╡ 2e8472f6-01b4-11eb-1f49-b35e9ee7c70b
model = DeBoerModel(TStar = 5000)

# ╔═╡ 373eadc4-01b6-11eb-0ef3-73d0e4518a24
md"""
N:
$(@bind N NumberField(1:1000, default = 100))
"""

# ╔═╡ 232d91d8-01b6-11eb-08d6-3343be95b77f
S,D,I,T,R = predict!(model, N, burnIn = 200)

# ╔═╡ 579f32aa-01b6-11eb-3a79-8f2207aeb03c
begin
	s = plot(S, lab = "", ylabel = "mmHg")
	plot!(D, lab = "", title = "BP", ylabel = "mmHg")
	t = plot(T, lab = "", title = "T", ylabel = "")
	i = plot(I, lab = "", title = "I", ylabel = "ms")
	r = plot(R, lab = "", title = "R", ylabel = "a.u.", xlab = "beats")
	plot(s,t,i,r, layout = (4,1), size = (700,800))
end

# ╔═╡ 5b0037b0-01b8-11eb-3ac0-27ba57005e81
begin
	fs = 1.25
	psdplot(S, lab = "S", fs = fs)
	psdplot!(D, lab = "D", fs = fs)
	psdplot!(R, lab = "R", fs = fs)
	psdplot!(I, lab = "I", fs = fs)
end

# ╔═╡ 0c881790-27f5-11eb-1786-cb457238f1ae
begin
	y = [S D T R]
	patPDC, patSpec, patCoh = pdc(y, metric = "diag")
	pdcplot(patPDC, patSpec, cnames = ["S" "D" "T" "R"])

end

# ╔═╡ ac3a7eac-6bba-11eb-37ec-49854b22053f
md"
# Simulating a phenylephrine injection

"

# ╔═╡ b312bc92-6bbf-11eb-131e-c946cf3eb23f
begin
	xbp = 0:30
	yms = 9 .* xbp
	scatter(S2 .- minimum(S2), I2 .- minimum(I2), lab = "")
	plot!(xbp, yms, color = :black, lab = "BRS 9 ms/mmHg", ylab = "ms", xlab = "mmHg")
end

# ╔═╡ 2dba5f40-6c47-11eb-3b29-47da9597609f
begin
	modelf = DeBoerModel()
	f, gain, phase, A = BRSf(modelf)
	plot(f, gain, lab = "", xlab = "f [c/b]", ylab = "BRS [ms/mmHg]", ylims = [0,18])
	res = hline!([mean(gain)], color  =  :black, lab = "mean BRS: $(round(mean(gain), digits = 2)) ms/mmHg")
end


# ╔═╡ Cell order:
# ╠═5c3096b0-01b6-11eb-2c01-5ba0aac27b70
# ╠═16771e7c-6c42-11eb-2b63-6f23d0ffa8d7
# ╠═a96c0a2a-01b3-11eb-3975-5953fd76dd7b
# ╠═0c720f5f-a8be-4d8e-a803-461118c3b2be
# ╠═2e8472f6-01b4-11eb-1f49-b35e9ee7c70b
# ╟─373eadc4-01b6-11eb-0ef3-73d0e4518a24
# ╠═232d91d8-01b6-11eb-08d6-3343be95b77f
# ╟─579f32aa-01b6-11eb-3a79-8f2207aeb03c
# ╟─5b0037b0-01b8-11eb-3ac0-27ba57005e81
# ╠═02590918-6c43-11eb-3f54-ab5c83dd1171
# ╠═0c881790-27f5-11eb-1786-cb457238f1ae
# ╟─ac3a7eac-6bba-11eb-37ec-49854b22053f
# ╠═53b158a2-6bbb-11eb-0fa7-87390c33cb64
# ╠═b312bc92-6bbf-11eb-131e-c946cf3eb23f
# ╠═2dba5f40-6c47-11eb-3b29-47da9597609f
