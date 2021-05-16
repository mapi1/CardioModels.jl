import ControlSystemIdentification.arxar
import ControlSystemIdentification.arx
import ControlSystemIdentification.InputOutputData
import ControlSystemIdentification.AbstractIdData
import ControlSystemIdentification.time1
import ControlSystemIdentification.parameter_covariance
import ControlSystemIdentification.ls
import ControlSystemIdentification.params2poly2
import ControlSystemIdentification.parameter_covariance
import ControlSystems.roots

function arx(d::AbstractIdData, na, nb; inputdelay = ones(Int, size(nb)), λ = 0, estimator = \, stochastic = false)
    y, u, h = time1(output(d)), time1(input(d)), sampletime(d)
    @assert size(y, 2) == 1 "arx only supports single output."
    # all(nb .<= na) || throw(DomainError(nb,"nb must be <= na"))
    na >= 0 || throw(ArgumentError("na must be positive"))
    size(nb) == size(inputdelay) || throw(ArgumentError("inputdelay has to have the same structure as nb"))
    y_train, A = getARXregressor(vec(y), u, na, nb, inputdelay = inputdelay)
    w = ls(A, y_train, λ, estimator)
    a, b = params2poly2(w, na, nb, inputdelay = inputdelay)
    model = tf(b, a, h)
    if stochastic
        local Σ
        try
            Σ = parameter_covariance(y_train, A, w, λ)
        catch e
            println(e)
            return minreal(model)
        end
        return TransferFunction(Particles, model, Σ)
    end

    return model
end

"""
    Σ = parameter_covariance(y_train, A, w, λ=0)
"""
function parameter_covariance(y_train, A, w, λ = 0)
    σ² = var(y_train .- A * w)
    iATA = if λ == 0
        inv(A'A)
    else
        ATA = A'A
        ATAλ = ATA + λ * I
        ATAλ \ ATA / ATAλ
    end
    iATA = (iATA + iATA') / 2
    Σ = σ² * iATA + sqrt(eps()) * Matrix(LinearAlgebra.I, size(iATA))
end

# function getφ(y, u, na, nb; N = nothing)
#     if N === nothing
#         N = length(y)
#     end
#     m = max(na, nb)+1
#     φy = toeplitz(y[m:N], y[m:-1:m-na])
#     φu = toeplitz(u[m:N], u[m:-1:m-nb])
#     φ = [-φy φu]
#     return φ
# end 
# +
# function parameter_covariance2(d, G, H, na, nb, nd)
#     nax = na + nd
#     nbx = nb + nd
#     y, u, h = time1(output(d)), time1(input(d)), sampletime(d)
#     φ = getφ(y, u, nax, nbx)
#     N = 1000
#     us = lsim(1/H, randn(N), 1:N)[1][:]
#     ys = lsim(G, randn(N), 1:N)[1][:]
#     inv(A'A)
# end

function arxar(d::InputOutputData, na::Int, nb::Union{Int, AbstractVector{Int}}, nd::Int;
    H::Union{TransferFunction, Nothing} = nothing,
    inputdelay      = ones(Int, size(nb)),
    δmin::Real      = 0.001,
    iterations::Int = 10,
    estimator       = \,
    verbose::Bool   = false,
    λ::Real         = 0,
    stochastic::Bool= false
)
    # Input Checking
    na >= 0 || throw(ArgumentError("na($na) must be positive"))
    all(nb .>= 0 )|| throw(ArgumentError("nb($nb) must be positive"))
    nd >= 1 || throw(ArgumentError("nd($nd) must be positive"))
    iterations >= 1 || throw(DomainError("iterations($iterations) must be >0"))
	δmin > 0 || throw(ArgumentError("δmin($δmin) must be positive"))
    ninputs(d) == length(nb) || throw(ArgumentError("Length of nb ($(length(nb))) must equal number of input signals ($(ninputs(d)))"))

    # 1. initialize H, GF and v to bring them to scope
    if H === nothing
        H = tf([0], [1], 1)
        iter = 0 # initialization
    else
        iter = 1 # initial noisemodel is known, H is present
    end
    GF = tf([0], [1], 1)
    dF, dH, v = 0, 0, 0
    # Iterate
    eOld    = 0
    δ       = δmin + 1
    timeVec = 1:length(d.y)
    sim(G,u) = lsim(G, u, timeVec)[1][:]
    while iter <= iterations && δ >= δmin
        # Filter input/output according to errormodel H, after initialization
        if iter > 0
            yF = sim(1/H, output(d)')
            if ninputs(d) == 1
                uF = sim(1/H, input(d)')
            else
                uF = fill(1.0, size(d.u)) 
                for i in 1:ninputs(d)
                    uF[i, :] = sim(1/H, d.u[i,:])
                end
            end
        else
            # pass unfiltered signal for initialization
            yF = copy(output(d))
            uF = copy(input(d))
        end
		dF = iddata(yF, uF, d.Ts)

		# 2. fit arx model
        GF = arx(dF, na, nb, estimator = estimator, λ = λ, inputdelay = inputdelay)

        # 3. Evaluate residuals
        v = residuals(GF, d)

        # 4. check if converged
        e = var(v)
        if eOld != 0
            δ = abs((eOld - e) / e)
        end
        verbose && println("iter: $iter, δ: $δ, e: $e")
        eOld = e

        # 5. estimate new noise model from residuals
        dH = iddata(v, d.Ts)
        H = ar(dH, nd, estimator = estimator)
        
        iter += 1
    end
    # assumption: GF and H are uncorrelated
    GF = arx(dF, na, nb, estimator = estimator, λ = λ, inputdelay = inputdelay, stochastic = stochastic)
    # total residuals e
    e = lsim(1/H, v, 1:length(v))[1][:]
    H = ar(dH, nd, estimator = estimator, stochastic = stochastic)
    return (G = GF, H = H, e = e)
end


# # Test data
# file = joinpath(pwd(), "data/paced_respiration/paced_respiation/data_PT/eroglu_PT_struct.mat")
# paced = import_paced(file)
# unsafe_comparisons(true)

# n = 10
# interval = paced.intervals[2]
# plot(interval)
# s = detrend(interval.SBP, p = 0)
# t = detrend(interval.RR, p = 0)
# r = detrend(interval.resp, p = 0)
# ampDead = fitBaselli(s, t, r, maxiter = 10,  na = n, nb = n, nc = n, stochastic = true)
# postDead = postprocess(ampDead, s, t, r)

# s, t, r = predict3(ampDead, 100)
# ribbonplot(s)
# plot(r)
# d = iddata(s, [t r], 1)
# G, H, e = arxar(d, n, [n, n], n)
# G2, H2, e = arxar(d, n, [n, n], n, stochastic = false)
# bodeplot(G)
# bodeplot!(G2)
# bodeplot(G2)
# minreal(G2)

# Gtest = G2[1,1]
# numvec(Gtest)
# ControlSystems.roots(numvec(Gtest))
# pole(Gtest)

# Gg = zpk([1±0.5], [2±0.2, 3±0.2], 1, 1)
# bodeplot(ss(Gg))
# pole(Gg)
# roots(Gg.num)
# Gtest.matrix[1]

# Gtest = tf([-0.00473 ± 0.0013, 0.000797 ± 0.0013, -0.00562 ± 0.0012], [ 1.0, -1.01 ± 0.066, 0.589 ± 0.08, -0.19 ± 0.059], 1)
# bodeplot(Gtest)

#### Test Example from Ljung ####
# N = 500
# A = tf([1, -0.5, 0.06], [1, 0, 0], 1)
# B = tf([1, 0.7], [1, 0,0], 1)
# G = minreal(B / A)
# C = tf([1, 0.95], [1, 0], 1)
# H = minreal(1 / (C * A))

# u = rand(Normal(0, 1), N)
# e = rand(Normal(0, (0.001)), N)
# sim(G, u) = lsim(G, u, 1:N)[1][:]
# y = sim(G, u)
# v = sim(H, e)
# yv = y.+ v# .+ rand(Normal(0, 1), N)
# d = iddata(yv, u, 1)
# Gest, Hest, resid = arxar(d, 2,2,1, stochastic = true)
