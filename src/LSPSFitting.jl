module LSPSFitting
using Base
using Formatting
using Dierckx
using LsqFit
using StatsBase
import DataFrames: DataFrame, describe, select, Not
using Statistics
using Printf
# ENV["MPLBACKEND"] = "MacOSX"

export fit_trace, fit_direct

function fit_trace(x, y; n::Int=3)
    efunc(t, p) = p[1] .* ((1.0 .- exp.(-t ./ p[2])) .^ n) .* exp.(-t ./ p[3])
    p0 = [minimum(y), 3.0, 15.0,]
    lb = [-1000.0, 0.2, 1.]
    ub = [0.0, 25., 100.]
    efit = curve_fit(efunc, x, y, p0) # , maxIter=10000)
    y2 = efunc(x, coef(efit))
    return y2, efit
end

function rolling_mean3(arr, n)
    so_far = sum(arr[1:n])
    out = zero(arr[n:end])
    out[1] = so_far
    for (i, (start, stop)) in enumerate(zip(arr, arr[n+1:end]))
        so_far += stop - start
        out[i+1] = so_far / n
    end
    return out
end

"""
    fit_direct(x, y, mindy)
    
Fit the currents in y with a simple equation

If there is a large derivative (larger than mindy),
do the fit anyway, but report with a flag for the trace

returns the flag (true if derivative is too big)
the fit waveform, and the fit report (fit.coef will have the fit values)

"""
function fit_direct(x, y; n::Int=1, mindy::Float64=0.75e-2, blank_window = [-1.0, 3.0])
    x = x .* 1e3  # msec
    y = y .* 1e12 # pA
    y_spl = Spline1D(x, y)
    dy = derivative(y_spl, x)
    yf = y
    # println("mindy: ", minimum(dy))
    minflag = false
    
    okpts = findall(abs.(rolling_mean3(dy, 7)) .< mindy)
    # println("mindy: ", mindy)
    # println("abs dy: ", abs.(dy))
    # println("okpts: ", okpts)
    y2, efit = fit_trace(x[okpts], y[okpts], n=n)
    
    # if minimum(dy) < mindy  # was -1e-7 for data in A
    #     minflag = true   # means peak value of derivative is too large - fit, but flag
    #     dt = mean(diff(x))
    #     pkmin = argmin(dy)
    #     pks = Int64(floor(pkmin + blank_window[1] / dt))
    #     pke = Int64(floor(pkmin + blank_window[2] / dt))
    #     if pks < 1
    #         pks = 1
    #     end
    #     if pke > size(y)[1]
    #         pke = size(y)[1]
    #     end
    #     yf = vcat(y[1:pks], y[pke:end])
    #     xf = vcat(x[1:pks], x[pke:end])
    #     y2, efit = fit_trace(xf, yf, n=n) # fit with "event" removed
    # else
    #     y2, efit = fit_trace(x, y, n=n) # no event to remove, so fit all
    # end
    efunc(t, p) = p[1] .* ((1.0 .- exp.(-t ./ p[2])) .^ n) .* exp.(-t ./ p[3])
    y0 = efunc(x, coef(efit))  # now compute comoplete y with fit parameters
    return minflag, y0*1e-12, efit  # note rescale... 
end

end
