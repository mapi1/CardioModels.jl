"""
psdplot(signal; fs = 1, maxf::Real = -1, method = "welch", order = 20)

Conviniently plots the power spectral density using one of two methods currently implemented in DSP.jl or a custom method using an AR model.
Futher plot parameters can be passed as keyword arguments.

# Args:

* 'signal::Vector': Data Vector containing te signal
* 'fs': Sample rate in Hz
* 'maxf': maximum frequency that shal be displayed
* 'method': Method for esitmation, either 'welch', 'ar' or 'fft'

# Examples

```jldoctest
julia> psd(signal)
```
"""
@userplot PSDPlot
# psdplot(signal, fs; maxf = nothing, estimator = welch_pgram, kwargs...)
@recipe function f(p::PSDPlot; fs = 1, maxf::Real = -1, method = "welch", order = 20)
    y = (p.args[1]) .- mean(p.args[1]) # detrend/ remove dc component
    # Methods for estimation
    if method == "ar"
        psd = pyulear(vec(y), order, fs = fs)
    elseif method == "fft"
        psd = periodogram(vec(y), fs = fs)
        # psd.power[1] = psd.power[2]
    elseif method == "welch"
        psd = welch_pgram(vec(y), fs = fs)
    else
        throw(ArgumentError("Invalid method choosen. Supported are 'welch', 'fft' and 'yule'."))
    end
    
    # plot porps
    yaxis := :log
    title := "PSD"
    # legend := false
    xaxis := "f [Hz]"
    if maxf > 0
        xlims := (0, maxf)
    else
        xlims := (minimum(psd.freq), maximum(psd.freq))
    end
    @series begin
        psd.freq, psd.power
    end
end

# Get PSD from AR process with variance σ² (Burg 1967)
function S(f::Number, σ::Number, φ::Vector{<:Number})
    indices = 1:length(φ)
    sumXes = sum(φ .* exp.(-2im * pi * f * indices))
    return 2σ^2 / abs2(1 + (sumXes))
end

"""
pyulear(x::Vector{<:Number}, order::Int)

Matlab Äquivalent von pyulear. Bestimme das Spektrum parametrisch mittels Yule-Walker-Gleichungen

# Args
* 'x::Vector{<:Number}': Zeitreihe aus der das Spektrum bestimmt werden soll
* 'order::Int': Ordnung des zugrundeliegenden AR Modells

# Keyword Args
* 'fs::Real = 1' Sampling frequency to get the power at the respective frequencies
"""
function pyulear(x::Vector{<:Real}, order::Int; fs::Real = 1, Δf = 0.001)
    f = 0:Δf:0.5
    fReal = range(0, stop = fs/2, length = length(f))
    # fit ar model
    (a,g) = lpc(x, order)
    # estimate spectrum
    Sar(x) = S(x, sqrt(g), a)
    power = Sar.(f)
    return DSP.Periodograms.Periodogram(power, fReal)
end


"""
getSpectralComponent(signal::Vector{<:Real}, f_range::Tuple{<:Real, <:Real}; method::String = "ar", Δf::Real = 0.001, order::Int = 10)

Get the spectral power of a certain frequency range, in natural frequencies, using different methods ['ar', 'fft', 'welch'].

# Args:

* 'signal::Vector{<:Real}': The signal that is to be decomposed.
* 'f_range::Tuple{<:Real, <:Real}': The frequency range of interest in natural frequencies [0, 0.5].

# Kwargs

* 'method::String = "ar"': Method used for psd  estimation. Choose from 'ar', 'fft' and 'welch'. 
* 'Δf::Real = 0.001': Only relevant for method ar, definfes distance of the grid of the underlying integration
* 'order::Int = 10': Only relevant for method ar, defines the order of the ar model 

# Examples

```jldoctest
julia> getSpectralComponent(signal, (0.1, 0.4))
```
"""
function getSpectralComponent(signal::Vector{<:Real}, f_range::Tuple{<:Real, <:Real}; method::String = "ar", Δf::Real = 0.001, order::Int = 10)
    0.0 <= f_range[1] < f_range[2] <= 0.5 || throw(DomainError("Start of range ($(f_range[1])) has to be smaller than end ($(f_range[2])) And in range of 0:0.5.")) 
    if method == "ar"
        psd = pyulear(signal, order, Δf = Δf)
    elseif method == "fft"
        psd = periodogram(signal)
        Δf = psd.freq[2] - psd.freq[1]
    elseif method == "welch"
        psd = welch_pgram(signal)
        Δf = psd.freq[2] - psd.freq[1]
    else
        throw(DomainError("Valid methods are: 'ar', 'fft' & 'welch'."))
    end
    pow = sum(Δf .* psd.power[findall(x -> f_range[1] <= x <= f_range[2], psd.freq)])
    return pow
end

"""
getSpectralComponentAR(ar::Union{ControlSystems.TransferFunction, AbstractVector}, f_range::Tuple{<:Real, <:Real}; σ::Real = 1,  Δf::Real = 0.001)

Get the spectral power of a certain frequency range, in natural frequencies, using different methods ['ar', 'fft', 'welch'].

# Args:

* 'ar::Union{ControlSystems.TransferFunction, AbstractVector}': Fitted AR model that schal be decomposed .
* 'f_range::Tuple{<:Real, <:Real}': The frequency range of interest in natural frequencies [0, 0.5].

# Kwargs

* 'Δf::Real = 0.001': Only relevant for method ar, definfes distance of the grid of the underlying integration
* 'σ::Real = 1': sqrt(innovation variance / prediction error)

# Examples

```jldoctest
julia> getSpectralComponent(ar(iddata(signal), 10), (0.1, 0.4))
```
"""
function getSpectralComponentAR(ar::Union{ControlSystems.TransferFunction, AbstractVector}, f_range::Tuple{<:Real, <:Real}; σ::Real = 1,  Δf::Real = 0.001)
    0.0 <= f_range[1] < f_range[2] <= 0.5 || throw(DomainError("Start of range ($(f_range[1])) has to be smaller than end ($(f_range[2])) And in range of 0:0.5.")) 
    if typeof(ar) <: ControlSystems.TransferFunction
        arcoef = denvec(ar)[1]
        arcoef = arcoef[2:end]
    elseif typeof(ar) <: AbstractVector
        arcoef = ar
    end
    f = 0:Δf:0.5
    Sar(x) = S(x, σ, arcoef)
    power = Sar.(f)
    pow = sum(Δf .* power[findall(x -> f_range[1] <= x <= f_range[2], f)])
    return pow
end