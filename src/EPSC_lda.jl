module EPSC_LDA
import DataFrames: DataFrame, describe, select, Not
import StatsBase: countmap, cor, var
using Distributions
using MLJ
using Distances
using PyPlot
export epsc_lda
"""
Compute LDA and identify events based on training set. 
d is the output of LPSPAnalysis.jl for a given run.

"""
function unpack_events(d)
    # unpack d from evts and return parallel arrays
    n = size(d.events)[1]
    trace = Vector{Int64}(undef, n)
    eventno = Vector{Int64}(undef, n)
    amp = Vector{Float64}(undef, n)
    peakt = Vector{Float64}(undef, n)
    onsettime = Vector{Float64}(undef, n)
    lat = Vector{Float64}(undef, n)  # relative to nearest previous 
    dur = Vector{Float64}(undef, n)
    rt = Vector{Float64}(undef, n)
    ft = Vector{Float64}(undef, n)
    rfratio = Vector{Float64}(undef, n)
    class = Vector{String}(undef, n)
    newclass = Vector{String}(undef, n)
    for i = 1:n
        trace[i] = d.events[i].trace
        eventno[i] = d.events[i].eventno
        amp[i] = d.events[i].amplitude
        peakt[i] = d.events[i].peaktime
        onsettime[i] = d.events[i].onsettime
        lat[i] = d.events[i].latency
        dur[i] = d.events[i].duration
        rt[i] = d.events[i].risetime
        rfratio[i] = rt[i] / d.events[i].falltime
        class[i] = d.events[i].class
        newclass[i] = d.events[i].newclass
        
    end
    return n, trace, eventno, amp, onsettime, peakt, lat, dur, rt, ft, rfratio, class, newclass
end

function events_to_dataframe(d)
    # unpack d and reformat as a DataFrame
    n, trace, eventno, amp, onset_t, peakt, lat, dur, rt, ft, rfratio, class, newclass =
        unpack_events(d)
    df = DataFrame(
        trace = trace,
        eventno = eventno,
        amp = amp,
        onsettime = onset_t,
        peaktime = peakt,
        lat = lat,
        dur = dur,
        rt = rt,
        rfratio = rfratio,
        class = class,
        newclass = newclass
    )
    return n, df
end
function epsc_lda(df)

    # n, amp, lat, dur, rt, ft, rfratio, class = unpack_events(d)
    # n, df = events_to_dataframe(d)
    r3 = x -> round(x, sigdigits = 3)

    y = df.class
    X = select(df, Not([:class, :newclass]))
    n = length(y)
    y = coerce(y, OrderedFactor)
    classes(y[1])
    # figure(figsize=(8,6))
    # cm = countmap(y)
    # bar([1, 2], [cm["spontaneous"], cm["direct"], cm["evoked"]])
    # xticks([1, 2], ["spontaneous", "direct", "evoked"], fontsize=12)
    # yticks(fontsize=12)
    # ylabel("Number of occurences", fontsize=14)

    train = 1:Int64(floor(n * 0.5))
    test = last(train)+1:n
    X2 = X
    # """
 #    Logistic regression model ("not bad")
 #    """
 #    println("LOGISTIC REGRESSION")
 #    @load LogisticClassifier pkg = MLJLinearModels
 #    X2 = X # select(X) # , Not([:]))
 #    clf = machine(LogisticClassifier(), X2, y)
 #    fit!(clf)
 #    ŷ = MLJ.predict(clf, X2)
 #    cross_entropy(ŷ, y) |> mean |> r3
 #    ŷ2 = predict_mode(clf, X2)
 #    println("Misclassification rate: ", misclassification_rate(ŷ2, y) |> r3)
 #    println("Accuracy: ", accuracy(ŷ2, y) |> r3)
 #    cm = confusion_matrix(ŷ2, y)
 #    display("text/plain", cm)


    # """
    # LDA on same data set
    # """
    # @load BayesianLDA pkg = MultivariateStats
    #
    # println("LDA with Bayesian")
    # clf = machine(BayesianLDA(), X2, y)
    # fit!(clf, rows = train)
    # ŷ_lda = predict_mode(clf, rows = test)
    #
    # println("Accuracy: ", accuracy(ŷ_lda, y[test]) |> r3)
    # cm = confusion_matrix(ŷ_lda, y[test])
    # display("text/plain", cm)

    """
    or: 
    """
    @load LDA pkg = MultivariateStats

    println("LDA")
    clf = machine(LDA(dist = Euclidean()), X2, y)  # CosineDist is; Euclidean seems to be better for accuracy
    mach = fit!(clf, rows = train)
    ŷ = predict_mode(clf, rows = test)
    println(ŷ[1:10])
    println(y[test[1:10]])
    println("Accuracy: ", accuracy(ŷ, y[test]) |> r3)
    cm = confusion_matrix(ŷ, y[test])
    display("text/plain", cm)
    # could do a save:
    MLJ.save("event_parser.jlso", mach)
    # then de-serialize:
    mach2 = machine("event_parser.jlso")
    #and run predict using the new classifier
    ynew = predict(mach2, X2)
    # and then ynew\hat will have the estimators (classifiers) for the new data
    # pdf(ynew, levels(ynew))
    z = predict_mode(mach2, X2)
    df.newclass[:] = z
    # for i = 1:20
    #     println(df.class[i], "   ", df.newclass[i])
    # end
    return df
    
    # """
    # Bayesian QDA from SciKit Learn
    # """
    # @load BayesianQDA pkg=ScikitLearn
    # BayesianQDA(
    #     priors = nothing,
    #     reg_param = 0.0,
    #     store_covariance = false,
    #     tol = 0.0001)
    # clf = machine(BayesianQDA(), X2, y)
    # fit!(clf, rows=train)
    # ŷ = predict_mode(clf, rows=test)
    #
    # accuracy(ŷ, y[test]) |> r3
end
end
