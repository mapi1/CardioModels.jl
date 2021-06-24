# Implementation of Karemakers adaption of the DeBoer Model
"""
    KaremakerModel()

A simple structure that stores a Karemaker-model. 

The default initialization was taken from the original paper. When constructiong the model, single properties can be defined as keyword arguments. 

# Examples
```julia
julia> model = KaremakerModel()
KaremakerModel
  c: Float64 350.0
  magn_int: Float64 10.0
  r_int: Distributions.Normal{Float64}
  A_resp: Float64 60.0
  T_resp: Float64 3000.0
  a1: Float64 9.0
  a2: Float64 9.0
  b_I: Float64 0.125
  b_D: Float64 0.125
  b_P: Float64 0.125
  Dref: Float64 80.0
  t_runoff: Float64 2000.0
  offset: Float64 10.0
  starling: Float64 0.016
  magn_pulse: Float64 2.0
  A_BP: Float64 3.0
  w: Array{Float64}((6,)) [0.0, 0.01, 0.03, 0.03, 0.015, 0.005]
  symp_set: Float64 0.5
  symp_add: Array{Float64}((2,)) [-0.05, 0.2]
```
"""
@with_kw mutable struct KaremakerModel
    c::Float64 = 350 # ms
    magn_int::Float64 = 10 # ms
    # r_int::Distribution = TriangularDist(-1, 1)
    r_int::Distribution = Normal(0, 1)
    # Respiration
    A_resp::Float64 = 60 # ms
    T_resp::Float64 = 3000 # ms
    # Baro 
    a1::Float64 = 9 # mmHg/ms
    a2::Float64 = 9 # mmHg/ms
    # sympatic force in I 
    b_I::Float64 = 0.125 
    # sympatic force in D 
    b_D::Float64 = 0.125 
    # sympatic force in P 
    b_P::Float64 = 0.125 
    Dref::Float64 = 80 # mmHg
    t_runoff::Float64 = 2000 # ms
    offset::Float64 = 10 # mmHg
    starling::Float64 = 0.016 # mmHg/ms
    magn_pulse::Float64 = 2 # ms
    A_BP::Float64 = 3 # mmHg 
    w::Vector{Float64} = [0.000,0.010, 0.030, 0.030, 0.015, 0.005]
    # w::Vector{Float64} = [1,2,3,3,2,1]
    # w::Vector{Float64} = [1,1,1,1,1,1]+
    ## Experimental ##
    symp_set = 0.5
    # symp_add = [-0.0, 0.0]
    symp_add = [-0.05, 0.2]
end

"""
    predict(model::KaremakerModel, n::Int; burnIn::Int = 0) 

Predicts 'N' values of the cardiovascular variables for the respective model. With 'burnIn' a certain number af values can be dropped in the beginning.
    
# Examples
```julia
 julia> S, D, P, I, Symp, ρ = predict(model, 100)  
 (S = [...], D = [...], P = [...], I = [...], Symp = [...], ρ = [...])   
 ```
 """
function predict(model::KaremakerModel, n::Int; burnIn::Int = 50)
    t = Vector{Float64}(undef, n + burnIn)
    trm = Vector{Float64}(undef, n + burnIn)
    I = Vector{Float64}(undef, n + burnIn)
    D = Vector{Float64}(undef, n + burnIn)
    S = Vector{Float64}(undef, n + burnIn)
    P = Vector{Float64}(undef, n + burnIn)
    Symp = Vector{Float64}(undef, n + burnIn)
    SympD = Vector{Float64}(undef, n + burnIn)
    
    history = 6
    for i in 1:history
        t[i] = 0
        trm[i] = model.c + model.magn_int * rand(model.r_int) + model.A_resp * sin((2π * t[i]) / model.T_resp + 1.5π)
        I[i] = 1000
        D[i] = model.Dref 
        S[i] = 160
        P[i] = 40
        Symp[i] = 0
        SympD[i] = 0
        
    end
    
    for i in (history+1):(n+burnIn)
        #timebase
        t[i] = t[i-1] + I[i-1]
        # constant + noise + rspiration for I 
        trm[i] = model.c + model.magn_int * rand(model.r_int) + model.A_resp * sin((2π * t[i]) / model.T_resp + 1.5π)
        
        # Symp[i] = 2/π * atan((sum(model.w .* SympD[i .- (1:length(model.w))]))/ (sum(model.w))) 
        # if Symp[i] > model.symp_set
        #     Symp[i] += model.symp_add[1] * Symp[i] 
        # else
        #     Symp[i] += model.symp_add[2] * Symp[i] 
        # end

        # normalize weights or not?
        Symp[i] = 2/π * atan((sum(model.w .* SympD[i .- (1:length(model.w))]))/ sum(model.w))
        # Symp[i] = atan((sum(model.w .* SympD[i .- (1:length(model.w))]))) # * 2/π    
        if I[i-1] > 700
            I[i] = ((model.a1 * S[i-1] + trm[i]) / (1 + model.b_I * Symp[i]))
        else
            I[i] = ((model.a2 * S[i-2] + trm[i]) / (1 + model.b_I * Symp[i]))
        end
        D[i] = S[i-1] * exp(-I[i] / (model.t_runoff * (1 + model.b_D * Symp[i]))) 
        P[i] = (model.offset + model.starling * I[i]) * (1 + model.b_P * Symp[i]) + model.magn_pulse * rand(model.r_int) + model.A_BP * sin((2π * t[i]) / model.T_resp)
        S[i] = D[i] + P[i]
        if D[i] > model.Dref
            SympD[i] = 0 
        else
            SympD[i] = model.Dref - D[i] 
        end
    end
    ρ = model.A_resp .* sin.((2π .* t) ./ model.T_resp .+ 1.5π)
    return ( S = S[burnIn+1:end], D = D[burnIn+1:end], P = P[burnIn+1:end], I = I[burnIn+1:end], Symp = Symp[burnIn+1:end], ρ = ρ[burnIn+1:end])# , SympD = SympD[burnIn+1:end])
    # return (I = I[burnIn+1:end], S = S[burnIn+1:end], D = D[burnIn+1:end], Symp = Symp[burnIn+1:end])
end


# model = KaremakerModel()
# i, s, d, symp = predict(model, 1000, burnIn = 45)
# get_lfhf(i)
# mean(symp)
# pi = plot(i)
# plot!(pi, i)
# ps = plot(s)
# plot!(ps, s)
# pd = plot(d)
# plot!(pd, d)
# lfhfs = []
# symps = []
# rr = []
# for runoff in 1250:2850
#     model = KaremakerModel(t_runoff = runoff)
#     i, s, d, symp = predict(model, 1000, burnIn = 50)
#     pgram = welch_pgram(i)
#     inds_lf = findall(x -> 0.04 <= x <= 0.15, pgram.freq)
#     inds_hf = findall(x -> 0.15 <= x <= 0.4, pgram.freq)
#     pow_lf = sum(pgram.power[inds_lf])
#     pow_hf = sum(pgram.power[inds_hf])
#     lfhf = pow_lf / pow_hf
#     push!(lfhfs, lfhf)
#     push!(lfhfs, lfhf)
#     push!(symps, mean(symp))
#     push!(rr, mean(i))
# end

# plot(symps, rr)

function get_lfhf(i)
    fs = 1000 / (mean(i))
    i = i .- mean(i)
    pgram = welch_pgram(i, fs = fs)
    inds_lf = findall(x -> 0.04 <= x <= 0.15, pgram.freq)
    inds_hf = findall(x -> 0.15 <= x <= 0.5, pgram.freq)
    pow_lf = sum(pgram.power[inds_lf])
    pow_hf = sum(pgram.power[inds_hf])
    lfhf = pow_lf / pow_hf  
    return lfhf
end


"""
    function vagal_balance(Dref::Real, runoffs::Union{AbstractVector{Real}, StepRangeLen}; n::Int = 5000, model::Union{Nothing, KaremakerModel} = nothing)

Calculates the LF/HF ratio for all distolic runoff constants in 'runoffs' returns also some means for different signals.
"""
function vagal_balance(Dref::Real, runoffs::Union{AbstractVector{Real}, StepRangeLen}; n::Int = 5000, model::Union{Nothing, KaremakerModel} = nothing)
    rr = []
    symp_drive = []
    lfhf = []
    sb = []
    db = []
    for runoff in runoffs
        if model  === nothing
            model = KaremakerModel(Dref = Dref, t_runoff = runoff)
        else
            model.Dref = Dref
            model.t_runoff = runoff
        end
        s, d, _, i, symp, _ = predict(model, n) 
        # i, s, d, symp = predict(model, n)
        push!(rr, mean(i))
        push!(sb, mean(s))
        push!(db, mean(d))
        push!(lfhf, get_lfhf(i))
        push!(symp_drive, mean(symp))
    end
    return rr, symp_drive, lfhf, sb, db
end


function vagal_balance_respiration(Ts::Union{AbstractVector{Real}, StepRangeLen}, Dref::Real, runoffs::Union{AbstractVector{Real}, StepRangeLen}; n::Int = 5000, model::Union{Nothing, KaremakerModel} = nothing)
    rr = []
    symp_drive = []
    lfhf = []
    for T in Ts
        if model  === nothing
            model = KaremakerModel(T_resp = T)
        else
            model.T_resp = T
        end
        rr_T, symp_drive_T, lfhf_T = vagal_balance(Dref, runoffs, model = model)
        push!(rr, rr_T)
        push!(lfhf, lfhf_T)
        push!(symp_drive, symp_drive_T)
    end
    return rr, symp_drive, lfhf
end

function lfhf2D_respiration(Ts::Union{AbstractVector{Real}, StepRangeLen}, Dref::Real, runoffs::Union{AbstractVector{Real}, StepRangeLen}; n::Int = 5000, model::Union{Nothing, KaremakerModel} = nothing)
    rr = []
    symp_drive = []
    lf = []
    hf = []
    for T in Ts
        if model  === nothing
            model = KaremakerModel(T_resp = T)
        else
            model.T_resp = T
        end
        rr_T, symp_drive_T, lf_T, hf_T = lfhf2D(Dref, runoffs, model = model)
        push!(rr, rr_T)
            push!(lf, lf_T)
        push!(hf, hf_T)
        push!(symp_drive, symp_drive_T)
    end
    return rr, symp_drive, lf, hf
end
