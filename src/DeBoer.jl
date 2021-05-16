# Implementation of the DeBoer Model

# include(joinpath(pwd(), "Tools/util.jl"))
# include(joinpath(pwd(), "Tools/tf2spec.jl"))

@with_kw mutable struct DeBoerModel
    # Gains for RR from Seff
    a0::Real = 9 # vagus (fast) ms/mmHG
    ak::Vector{<:Real} = [0, 1, 2, 3, 2, 1] # sympathicus (slow) ms/mmHG
    c1::Real = 800 - 120a0 - sum(120 .* ak) # constant for mean RR interval
    # Gains for peripheral resistance
    bk::Vector{<:Real} = [0, 2, 4, 6, 4, 2]
    TStar::Real = 1425 + sum(120 .* bk) # Base aterial time constant
    # Baroreceptor Gain
    γ::Real = 0.016 # mmHg/ms
    c2::Real = 120 - 75 - 800γ
    #Windkessel
    c3::Real = 75 / 120exp(-800 / 1425)
    # Respiration
    A::Real = 3 # mmHg Amplitude of respirational influence
    fresp::Real = 0.3 # 18 breaths/min
    resp::Function = t -> A*sin(2pi * fresp * t / 1000) # Sinus Respiration
    # Noise (assumed as normally distributed with 0 mean)
    noiseI::Distribution = Normal(0, 25) # ms
    noiseS::Distribution = Normal(0, 2)  # mmHg
    # Function for effective BP Seff
    operationPoint::Real = 120
    F::Function = S -> operationPoint .+ 18atan.((S .- operationPoint) ./ 18)

    # model state
    hasState::Bool = false
    S::Union{Nothing, Vector{<:Real}} = nothing
    Seff::Union{Nothing, Vector{<:Real}} = nothing
    D::Union{Nothing, Vector{<:Real}} = nothing
    I::Union{Nothing, Vector{<:Real}} = nothing
    T::Union{Nothing, Vector{<:Real}} = nothing
    R::Union{Nothing, Vector{<:Real}} = nothing
    time::Real = 0
end

# Statefull predict function
function predict!(model::DeBoerModel, n::Int; burnIn::Int = 0)
    burnIn >= 0 || throw(DomainError(burnIn, "Burn in needs to >= 0"))
    n > 0 || throw(DomainError(n, "You need to predict at least one beat"))

    lena = length(model.ak) + 1
    lenb = length(model.bk)
    history = max(lena, lenb)
    if model.hasState && history <= length(model.S)
        # Initialization Phase - connect new arrays to state
        S = append!(model.S, Vector{Real}(undef, n + burnIn))
        Seff = append!(model.Seff, Vector{Real}(undef, n + burnIn))
        D = append!(model.D, Vector{Real}(undef, n + burnIn))
        I = append!(model.I, Vector{Real}(undef, n + burnIn))
        T = append!(model.T, Vector{Real}(undef, n + burnIn))
        R = append!(model.R, Vector{Real}(undef, n + burnIn))
        time = model.time
        burnIn += history
    else
        # Default Initialization
        if model.hasState
            @warn "State was too short (length of $(length(model.S)) but history needed is $history), switching to default initialization"
        end
        burnIn += history
        S = Vector{Real}(undef, n + burnIn)
        Seff = Vector{Real}(undef, n + burnIn)
        D = Vector{Real}(undef, n + burnIn)
        I = Vector{Real}(undef, n + burnIn)
        T = Vector{Real}(undef, n + burnIn)
        R = Vector{Real}(undef, n + burnIn)
        time = 0
        # fill in starting data
        for i in 1:history
            S[i] = model.operationPoint
            Seff[i] = model.F(model.operationPoint)
            D[i] = 75
            I[i] = 800
            T[i] = 1425
        end
    end

    # generate future data
    for i in history+1:n+burnIn
        time += I[i-1]
        D[i] = model.c3 * S[i-1] * exp(-I[i-1]/T[i-1])
        R[i] = model.resp(time)
        S[i] = D[i] + model.γ * I[i-1] + model.c2 + R[i] + rand(model.noiseS)
        Seff[i] = model.F(S[i])
        I[i] = model.c1 + model.a0 * Seff[i] + sum(model.ak .* reverse(Seff[i-(lena-1):i-1])) + rand(model.noiseI)
        T[i] = model.TStar - sum(model.bk .* reverse(Seff[i-lenb:i-1]))
    end

    # update / shorten  state
    model.hasState = true
    model.S = S[end-history+1:end]
    model.Seff = Seff[end-history+1:end]
    model.D = D[end-history+1:end]
    model.I = I[end-history+1:end]
    model.T = T[end-history+1:end]
    model.R = R[end-history+1:end]
    model.time = time

    return (S = S[burnIn+1:end], D = D[burnIn+1:end], I = I[burnIn+1:end], T = T[burnIn+1:end], R = R[burnIn+1:end])
end

# simulates a phenylephrine injection by increaing the peripheral resistance over a defined time 
function phenylephrine(model::DeBoerModel, TStar_increase::Real, n_increase::Int)
    n_increase > 0 || throw(DomainError(n_increase, "peripheral resistance needs to invrease over a minimum of one beat"))
    TStar_step = TStar_increase / n_increase
    
    S = Vector{Float64}(undef, n_increase)
    D = Vector{Float64}(undef, n_increase)
    I = Vector{Float64}(undef, n_increase)
    T = Vector{Float64}(undef, n_increase)
    R = Vector{Float64}(undef, n_increase)
    for i in 1:n_increase
        model.TStar += TStar_step
        res = predict!(model, 1)
        S[i] = res[:S][1]
        D[i] = res[:D][1]
        I[i] = res[:I][1]
        T[i] = res[:T][1]
        R[i] = res[:R][1] 
    end
    return (S = S, D = D, I = I, T = T, R = R)
end

function BRSf(model::DeBoerModel)
    a = [model.a0; model.ak]
    A = impRes2tf(a)
    w = range(0, stop = π, length = 256)
    mag, phase, w = bode(A, w)
    f = w ./ 2π
    return f, vec(mag), vec(phase), A
end

## BRS(f)
# model = DeBoerModel()
# f, gain, phase, A = BRSf(model)
# plot(f, gain, lab = "", xlab = "f [c/b]", ylab = "BRS [ms/mmHg]", ylims = [0,18])
# hline!([9], color  =  :black, lab = "vagal BRS: 9 ms/mmHg")
##

### f(S)
# model = DeBoerModel()
# S = 60:180
# Seff = model.F.(S)
# plot(S, Seff, ylab = "Seff", xlab = "S", lab = "", color = :black, size = (200, 200))
# savefig("./results/seff.svg")
###

### sigmoid for BRS
# sig(t) = 1 / (1 + exp(-t))
# t = -5:0.01:5
# plot(t, sig.(t), lab = "", grid = :off, color = :black, linewidth = 2)
# savefig("./results/brs.svg")
# x = 1
###
# x = 1
# s,d,i,t = predict(model, 100)
# plot([s, i])