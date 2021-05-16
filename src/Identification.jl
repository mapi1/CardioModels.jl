include("stochastic_arxar.jl")

# fit model from data
function fitBaselli(S::Vector{<:Real}, I::Vector{<:Real}, ρ::Vector{<:Real};
     na::Int = 10,
     nb::Int = 10,
     nd::Int = 10,
     maxiter::Int = 10,
     δmin::Real = 0.001,
     verbose::Bool = false,
     λ::Union{Vector{<:Real}, Real} = 0,
     stochastic::Bool = false)
    
     # Check arguments
    length(S) == length(I) == length(ρ) || throw(ArgumentError("S, I, ρ must be of the same length"))
    na >= 1 || throw(ArgumentError("na ($na) must be positive"))
    nb >= 1 || throw(ArgumentError("nb ($nb) must be positive"))
    nd >= 1 || throw(ArgumentError("nd ($nd) must be positive"))
    maxiter >= 1 || throw(ArgumentError("maxiter($maxiter) must be positive"))
    δmin >= 0 || throw(ArgumentError("δmin($δmin) must be positive"))
    
    # Prepare inputs
    N = length(S)
    T = 1

    us = [I ρ]
    ds = iddata(S, us, T)

    ut = [S ρ]
    dt = iddata(I, ut, T)
    
    # Estimate Transfer Functions using generalized least squares, either the maximum likelihood version or repeated arxar 
    # ARXXAR
    # Gs, Ds, errs, residual_s = gls(ds, na, [1, nb], nd, δmin = δmin, maxiter = maxiter, errorAndResidual = true, verbose = verbose, estimator = estimator,  λ = λ)
    Gs, Ds, residual_s = arxar(ds, na, [1, nb], nd, δmin = δmin, iterations = maxiter, verbose = verbose, estimator = estimator,  λ = λ, stochastic = stochastic)

    # XXAR - na = 0 
    # Gt, Dt, errt, residual_t = gls(dt, 0, [nb, nb], nd, direct = [true, false], δmin = δmin, maxiter = maxiter, errorAndResidual = true, verbose = verbose, estimator = estimator,  λ = λ)
    Gt, Dt, residual_t = arxar(dt, 0, [nb, nb], nd, inputdelay = [0, 1], δmin = δmin, iterations = maxiter, verbose = verbose, estimator = estimator,  λ = λ, stochastic = stochastic)

    # Extract model
    # G = G_SI / (1 - G_SS) with a priori knowledge gSI(1) = g(1) gives G_SI and G_SS
    gSI = vec(impulse(G,1)[1]')
    G_SI = impRes2tf(gSI)
    G_SS = ControlSystems.robust_minreal(-(G_SI / G) + 1)
    R1 = Gs[1,2]
    R_Sρ = ControlSystems.robust_minreal(R1 * (1 - G_SS))
    MS = ControlSystems.robust_minreal(Ds * (1 - G_SS))
    
    G_IS = Gt[1,1]
    R_Iρ = Gt[1,2]
    MI = Dt
    
    # Fit AR model representing respiration
    arCoeffs, errr = lpc(ρ, nd)
    Mρ = tf(1, ([1; arCoeffs]), 1)
    # TODO use ar(;stochastic = stochastic)

    # Check residuals for witheness(normal distribution) via 5% anderson darling test
    pS_anderson = pvalue(OneSampleADTest(residual_s, Normal(mean(residual_s), std(residual_s))))
    pI_anderson = pvalue(OneSampleADTest(residual_t, Normal(mean(residual_t), std(residual_t))))
    
    # res = BaselliModel(G_SI = G_SI, G_IS = G_IS, G_SS = G_SS, MS = MS, MI = MI, Mρ = Mρ, R_Sρ = R_Sρ, R_Iρ = R_Iρ, wS = Normal(0, sqrt(errs)), wI = Normal(0, sqrt(errt)), wρ = Normal(0, sqrt(errr)), n = na, N = N, pS_anderson = pS_anderson, pI_anderson = pI_anderson, res_S = residual_s, res_I = residual_t, Gs = Gs, Gt = Gt, Ds = Ds,  Dt = Dt)
    res = BaselliModel(G_SI = G_SI, G_IS = G_IS, G_SS = G_SS, MS = MS, MI = MI, Mρ = Mρ, R_Sρ = R_Sρ, R_Iρ = R_Iρ, wS = Normal(0, std(residual_s)), wI = Normal(0, std(residual_t)), wρ = Normal(0, sqrt(errr)), n = na, N = N, pS_anderson = pS_anderson, pI_anderson = pI_anderson, res_S = residual_s, res_I = residual_t, Gs = Gs, Gt = Gt, Ds = Ds,  Dt = Dt)
    return res#, (Gs, Gt, Ds,  Dt)
end

@with_kw mutable struct rstData
    S::Vector{<:Real} = [1]
    I::Vector{<:Real} = [1]
    ρ::Vector{<:Real} = [1]
end

# extract physiologically relevant parameters from model
function postprocess(M::BaselliModel, data::rstData)
    # BRS
    # ramp = 1:140
    # rampRes = lsim(M.G_IS, ramp, ramp)[1][:]
    # α = theilSenEstimator(ramp, rampRes)[1]
    gIS = impulse(M.G_IS)[1]
    α = sum(gIS) # "BRS"
    a0 = gIS[1] # vagal part of BR
    a1 = gIS[2] # vagal part of BR
    
    # gSI stroke volume effect (+) or run off effect (-)
    gSI = vec(impulse(M.G_SI,1)[1]')[2]
    
    ## AP regulation
    # gSS represents atrial resistance and compliance (Windkessel)
    gSS = vec(impulse(M.G_SS,1)[1]')[2]
    # Gss(LF), resonance frequency of closed G_SS loop -> index of contribution of AP regulation LF oscilation
    Gss = 1 / (1 - M.G_SS)
    mag, phase, w = bode(Gss, range(0.0, stop = π, step = 0.0100))
    GssLF = maximum(mag)
    
    ## Respiratory effects - Share of HF power due to respiration / HF unwell defined for lof f_resp
    # Effect of resipration on I through different pathways
    # TF = minreal((1-M.G_SS)*M.R_Iρ + M.G_IS * M.R_Sρ) * Mρ / (1 - M.G_SS - M.G_SI * M.G_IS)
    # 1) R_Sρ is set to 0
    TF1 = ControlSystems.robust_minreal((1-M.G_SS)*M.R_Iρ + M.G_IS * 0) * M.Mρ / (1 - M.G_SS - M.G_SI * M.G_IS)
    # 2) R_Iρ is set to 0
    TF2 = ControlSystems.robust_minreal((1-M.G_SS)*0 + M.G_IS * M.R_Sρ) * M.Mρ / (1 - M.G_SS - M.G_SI * M.G_IS)  
    wHF = range(0.15*2π, stop = 0.4*2π, step = 0.001)
    # fHF = wHF / 2π
    mag1, phase, w = bode(TF1, wHF)
    mag2, phase, w = bode(TF2, wHF)
    # Modulus squared?
    Pt1 =  mag1[:] .* (M.wρ.σ ^2)
    Pt2 =  mag2[:] .* (M.wρ.σ ^2)
    Θt = log10(sum(Pt1) / sum(Pt2))
    
    # Effect of respiration on S through different pathways
    # Tf = (M.G_SI * M.R_Iρ + M.R_Iρ) * Mρ / (1 - M.G_SS - M.G_SI * M.G_IS)
    TF1 = (M.G_SI * M.R_Iρ + 0) * M.Mρ / (1 - M.G_SS - M.G_SI * M.G_IS)
    TF2 = (M.G_SI * 0 + M.R_Sρ) * M.Mρ / (1 - M.G_SS - M.G_SI * M.G_IS)
    mag1, phase, w = bode(TF1, wHF)
    mag2, phase, w = bode(TF2, wHF)
    # Modulus squared?
    Pt1 =  mag1[:] .* (M.wρ.σ ^2)
    Pt2 =  mag2[:] .* (M.wρ.σ ^2)
    Θs = log10(sum(Pt1) / sum(Pt2))
    
    # LF (0.04:0.15)sources in residuals
    LF = (0.04, 0.15)
    sLFpow = getSpectralComponentAR(M.MS, LF, σ = M.wS.σ)
    PusLF = sLFpow / var(detrend(data.S)) 
    tLFpow = getSpectralComponentAR(M.MS, LF, σ = M.wS.σ)
    PutLF = tLFpow / var(detrend(data.I))       

    return (α = α, a0 = a0, a1 = a1,gSI = gSI, gSS = gSS, GssLF = GssLF,  Θt =  Θt, Θs = Θs, PusLF = PusLF, PutLF = PutLF)
end

# wrapper without rstData
function postprocess(M::BaselliModel, S::Vector{<:Real}, I::Vector{<:Real}, ρ::Vector{<:Real})
    data = rstData()
    data.ρ = ρ
    data.I = I
    data.S = S
    return postprocess(M, data)
end

function minimal_postprocess(M::BaselliModel, S::Vector{<:Real}, I::Vector{<:Real}, ρ::Vector{<:Real})
    gIS = impulse(M.G_IS)[1]
    α = sum(gIS) # "BRS"
    return α
end
# Simple Parameter of Goodness for fitted model (R² & R²adjusted)
function getGoodness(M::BaselliModel, data::rstData)
    # R²
    vars = var(detrend(data.S))
    vart = var(detrend(data.I))
    ps = ((vars - (M.wS.σ ^2)) / vars) 
    pt = ((vart - (M.wI.σ ^2)) / vart)
    # R²adjusted
    ps_adjusted = 1 - ((M.N - 1) / (M.N - M.n)) * (1 - ps)
    pt_adjusted = 1 - ((M.N - 1) / (M.N - M.n)) * (1 - pt)
    return (ps = ps, pt = pt, ps_adjusted = ps_adjusted, pt_adjusted = pt_adjusted)
end

# wrapper without rstData
function getGoodness(M::BaselliModel, S::Vector{<:Real}, I::Vector{<:Real}, ρ::Vector{<:Real})
    data = rstData()
    data.ρ = ρ
    data.I = I
    data.S = S
    return getGoodness(M, data)
end
