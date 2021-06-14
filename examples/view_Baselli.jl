### A Pluto.jl notebook ###
# v0.14.3

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

# ╔═╡ 6d960bb1-b4d7-49f3-8549-78786eabfc74
using PlutoUI, Plots, KardioUtils

# ╔═╡ 49ac8ea6-a5dd-11eb-1952-19922f814f4d
include("./Baselli.jl")

# ╔═╡ d494e27e-7c67-4abc-b491-0dd09e58f169
md"
# Baselli Model
Modelling cardiorespiratory varibilities
"

# ╔═╡ e7684c53-fc26-4658-9825-d622c442f394
md"""
Simulation: 
$(@bind sim NumberField(1:4, default = 1))
"""

# ╔═╡ 8ad5bc71-dea2-43cc-b745-9c47ff6d7dff
model = getModel2(sim);

# ╔═╡ 8c387e69-c2c6-4e93-897a-bdd47c65764d
plot(model, lab = "")

# ╔═╡ 8f948c31-0561-42e6-a0eb-48891126d9aa
S, I, ρ = predict3(model,  512);

# ╔═╡ b4c65d6d-72a8-470d-96c3-3d0c6bd810aa
md"""
PSD: 
$(@bind ispsd CheckBox())
"""

# ╔═╡ 0993d11c-c05f-48e7-9c1d-c95bba3d7251
begin
	if ispsd
		pS = psdplot(S, lab = "", title = "")
		pI = psdplot(I, lab = "")
		pρ = psdplot(ρ, lab = "")
	else
		pS = plot(S, lab = "", ylab = "S [mmHg]")
		pI = plot(I, lab = "", ylab = "I [ms]")
		pρ = plot(ρ, lab = "", ylab = "ρ [a.u.]", xlab =  "beat #")
	end
	plot(pS, pI, pρ, layout = (3,1))
end

# ╔═╡ Cell order:
# ╠═6d960bb1-b4d7-49f3-8549-78786eabfc74
# ╠═49ac8ea6-a5dd-11eb-1952-19922f814f4d
# ╟─d494e27e-7c67-4abc-b491-0dd09e58f169
# ╟─e7684c53-fc26-4658-9825-d622c442f394
# ╠═8ad5bc71-dea2-43cc-b745-9c47ff6d7dff
# ╠═8c387e69-c2c6-4e93-897a-bdd47c65764d
# ╠═8f948c31-0561-42e6-a0eb-48891126d9aa
# ╟─b4c65d6d-72a8-470d-96c3-3d0c6bd810aa
# ╟─0993d11c-c05f-48e7-9c1d-c95bba3d7251
