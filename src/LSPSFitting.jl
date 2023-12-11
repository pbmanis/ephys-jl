module LSPSFitting
using AddPackage
using Base
using Formatting
using Dierckx
using LsqFit
using StatsBase
@add using NoiseRobustDifferentiation
import DataFrames: DataFrame, describe, select, Not
using Statistics
using Printf
# ENV["MPLBACKEND"] = "MacOSX"

export fit_trace, fit_direct

function direct_template(t, p)
    f = p[1] .* ((1.0 .- exp.(-(t .- p[4]) ./ p[2])) .^ n) .* exp.(-(t .- p[4]) ./ p[3])
    return f
end
"""
    fit_trace(x, y, n)
"""
function fit_trace(x, y; n::Int=2)
    # efunc(t, p) = p[1] .* ((1.0 .- exp.(-(t .- p[4]) ./ p[2])) .^ n) .* exp.(-(t .- p[4]) ./ p[3])
    p0 = [minimum(y), 3.0, 15.0, 0.]   # amplitude, rise time, decay time, delay
    lb = [-2000.0, 0.2, 1.0, 0.]
    ub = [0.0, 25.0, 100.0, 1.0]
    efit = curve_fit(direct_template, x, y, p0, lower=lb, upper=ub) # , maxIter=10000)
    y2 = direct_template(x, coef(efit))
    # print(coef(efit))
    return y2, efit
end

function rolling_mean3(arr, n)
    so_far = sum(arr[1:n])
    out = zero(arr[n:end])
    out[1] = so_far
    for (i, (start, stop)) in enumerate(zip(arr, arr[(n + 1):end]))
        so_far += stop - start
        out[i + 1] = so_far / n
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
function fit_direct(x, y; n::Int=1, mindy::Float64=0.75e-2, fitmode::String = "spline", blank_window=[-1.0, 3.0])
    x = x .* 1e3  # msec
    y = y .* 1e12 # pA
    # y_spl = Spline1D(x, y)
    # compute derivatives for subtraction
    dy = tvdiff(y, 200, 0.5; dx=x[2]-x[1], scale="small", ε = 1e-9)
    #dy = derivative(y_spl, x)

    # println("mindy: ", minimum(dy))
    minflag = false

    # find points where slope is less than some min value - these are "ok"
    okpts = findall(abs.(rolling_mean3(dy, 7)) .< mindy)
    println("npts: ", length(x), " n okpts: ", length(okpts))
    
    if fitmode === "spline" # smooth function with spline
        fit = Spline1D(x[okpts], y[okpts]) # refit
        template = evaluate(fit, x) # re-evaluate at all points

    else # or fit to a function
        template, fit = fit_trace(x[okpts], y[okpts]; n=n)
        # efunc(t, p) = p[1] .* ((1.0 .- exp.(-(t .- p[4]) ./ p[2])) .^ n) .* exp.(-(t .- p[4]) ./ p[3])
        template = direct_template(x, coef(fit))  # now compute complete y with fit parameters
    end
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

    # y0 = evaluate(spline_okpts, x) # re-evaluate at all points
    return minflag, template*1e-12, fit  # note rescale... 
end

end
