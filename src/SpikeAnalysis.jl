__precompile__()
module SpikeAnalysis
using LsqFit
using Statistics
using Printf

#=
This module provides functions for performing analysis of spikes.
1. spike detection: 3 methods:
	a. voltage threshold (with minimum width)
	b. slope and voltage-dependent (Hight and Kalluri)
	c. peak (default)
This module roughly duplicates the ephys.SpikeAnalysis class in python

NOTE: All times are in seconds.
      All voltages are in Volts.
	  All currents are in Amperes.

Data structures are used to simplify handling the results.
Voltagetrace holds the data for a single trace
Spike holds the analysis results for a single spike
Spikes holds an array of individual Spike data
SpikeDetection holds the parameters for the detection algorithm
=#

struct Spike
    indices::Tuple{Int64,Int64}
    trace::Int64 # trace of origin
    eventno::Int64 # spike number in trace
    latency::Float64 # latency in seconds
    amplitude::Float64  # peak kamplitude in mV
    onsettime::Float64
    peaktime::Float64
    risedvdt::Float64
    halfwidth::Float64
    falldvdt::Float64
    ahpdepth::Float64
    mode::String
    detector::String
end

struct Spikes
    label::String
    events::Vector{Spike}
end

mutable struct SpikeDetectParameters
    threshold::Float64 # threshold in V
    refractory::Float64# refractory period in sec
    min_halfwidth::Float64# minimum halfwidth in sec
    HK_dt2::Float64# specific parameters for Hight & Kalluri detector
    HK_C1::Float64 # V
    HK_C2::Float64 # V
    mode::String  # detection mode (schmitt, threshold or peak)
    detector::String  # detector method to use (threshold, argrelmax, Kalluri)

end
function SpikeDetectParameters(;
    threshold = 0.0,
    refractory = 0.0007,
    min_halfwidth = 0.0001,
    HK_dt2 = 0.00175,
    HK_C1 = -0.012,
    HK_C2 = 0.011,
    mode = "threshold",
    detector = "Kalluri",
)
    return (SpikeDetectParameters(
        threshold,
        refractory,
        min_halfwidth,
        HK_dt2,
        HK_C1,
        HK_C2,
        mode,
        detector,
    ))
end

function AnalyzeSpikes(
    tdat,
    vdat,
    idat;
    timewindow = [0.1, 0.6],
    pars = Union{SpikeDetectParameters,missing},
)
    # vtrace = VoltageTrace(tdat, vdat, idat, 0)
    # println(pars)
    if ismissing(pars)
        pars = SpikeDetectParameters()
    end
    if pars.detector == "Kalluri"
        spkt, spkv, spki = HightKalluri(tdat, vdat, idat, pars)
        # println(spkt, spki)
    end
    return spkt, spkv, spki
end


function HightKalluri(tdat, vdat, idat, pars)
    #=
    Find spikes using a box method:
    Voltage must be > threshold, and have slope values restricted to a range
    Units must be consistent: x, dt, d2 (s or ms)
    Unist must be consistent: y, thr, C1, C2 (V or mV)
    Note: probably works best with mV and ms, given the constants above.
    to C1, C2 and the width dt2
    From Hight and Kalluri, J Neurophysiol., 2016
    Returns an array of indices in x where spikes occur
    =#
    # parameters for the Hight and Kalluri detecctor
    dt2 = pars.HK_dt2  # s
    C1 = pars.HK_C1  # V
    C2 = pars.HK_C2  # V
    threshold = pars.threshold

    npts = size(vdat)
    spikes = zero(vdat)
    dt = tdat[2] - tdat[1] # sample rate
    iwid = round(Int64, dt2 / dt) # width of region to look for slope values
    for i = iwid:(length(vdat)-iwid)
        if vdat[i] > threshold # when voltage exceeds threshold, check shape
            if (vdat[i] > vdat[i-1]) & (vdat[i] > vdat[i+1])  # local peak
                if ((vdat[i+iwid] - vdat[i]) < C1) & ((vdat[i] - vdat[i-iwid]) > C2)
                    spikes[i] = 1.0
                end
            end
        end
    end
    spki = @view(spikes[2:end]) - @view(spikes[1:end-1]) .> 0
    insert!(spki, 1, 0) # offset for first element
    spike_indices = findall(x -> x == 1, spki)
    spkt = (spike_indices .- 1) .* dt
    spkv = vdat[spike_indices]
    return spkt, spkv, spike_indices
end

# function findspikes(vtrace, pars)
#     if pars.mode not in ["schmitt", "threshold", "peak"]:
#         raise ValueError(
#             'pylibrary.utility.findspikes: mode must be one of "schmitt", "threshold", "peak" : got %s'
#             % mode
#         )
#     if pars.detector not in ["threshold", "argrelmax", "Kalluri"]:
#         raise ValueError(
#             'pylibrary.utility.findspikes: mode must be one of "argrelmax", "threshold" "Kalluri": got %s'
#             % detector
#         )
#
#     if t1 is not None and t0 is not None:
#         xt = ma.masked_outside(x, t0, t1)
#         vma = ma.array(v, mask=ma.getmask(xt))
#         xt = ma.compressed(xt)  # convert back to usual numpy arrays then
#         vma = ma.compressed(vma)
#     else:
#         xt = np.array(x)
#         vma = np.array(vma)
#     # print('max x: ', np.max(xt))
#    #  print('dt: ', dt)
#
#     dv = np.diff(vma) / dt  # compute slope
#     dv2 = np.diff(dv) / dt  # and second derivative
#     st = np.array([])  # store spike times
#
#     if detector == "threshold":
#         spv = np.where(vma > thresh)[0].tolist()  # find points above threshold
#         sps = (
#             np.where(dv > 0.0)[0] + 1
#         ).tolist()  # find points where slope is positive
#         sp = list(
#             set(spv) & set(sps)
#         )  # intersection defines putative spike start times
#         # then go on to mode...
#
#     elif detector == "argrelmax":
#         #  spks = scipy.signal.find_peaks_cwt(vma[spv], np.arange(2, int(peakwidth/dt)), noise_perc=0.1)
#         order = int(refract / dt) + 1
#         stn = scipy.signal.find_peaks(vma, height=thresh, distance=order)[0]
#         # argrelmax seems to miss peaks occasionally
#         # spks = scipy.signal.argrelmax(vma, order=order)[0]
#         # stn = spks[np.where(vma[spks] >= thresh)[0]]
#         if len(stn) > 0:
#             stn2 = [stn[0]]
#         else:
#             stn2 = []
#         # filter peaks by checking that valleys between pairs
#         # are sufficiently deep. Note that this only checks
#         # BETWEEN spikes, so we need to do an additional
#         # check of the last "spike" separately
#         removed = []
#         t_forward = int(0.010 / dt)  # use 10 msec forward for drop
#         for i in range(len(stn) - 1):  # for all putative peaks
#             if i in removed:  # this can happen if event was removed in j loop
#                 continue
#             test_end = min([stn[i] + t_forward, stn[i + 1], vma.shape[0]])
#
#             if stn[i] == test_end:
#                 continue
#             elif (vma[stn[i]] - np.min(vma[stn[i] : test_end])) < mindip:
#                 if (
#                     i == 0
#                 ):  # special case: if first event fails, remove it from output list
#                     stn2 = []
#                 removed.append(i)
#                 continue
#             else:
#                 stn2.append(stn[i])
#         # handle "spikes" that do not repolarize and are the *last* spike
#         if len(stn2) > 1:
#             test_end = stn2[-1] + t_forward
#             minv = np.min(vma[stn2[-1] : test_end])
#             if (vma[stn2][-1] - minv) < mindip:
#                 removed.append(stn2[-1])
#                 stn2 = stn2[:-1]  # remove the last spike
#         stn2 = sorted(list(set(stn2)))
#         if debug:
#             print("stn: ", stn)
#             print(vma[stn])
#         # if len(stn2) > 0:  # good to test algorithm
#         #     import matplotlib.pyplot as mpl
#         #
#         #     f, ax = mpl.subplots(1,1)
#         #     ax.plot(xt, vma)
#         #     ax.plot(xt[stn2], vma[stn2], 'ro')
#         #     ax.plot(xt[stn], vma[stn], 'bx')
#         #     mpl.show()
#         xspk = x[[s + int(t0 / dt) for s in stn2]]
#         return self.clean_spiketimes(xspk, mindT=refract)  # done here.
#     else:
#         raise ValueError("Utility:findspikes: invalid detector")
#
#     sp.sort()  # make sure all detected events are in order (sets is unordered)
#
#     spl = sp
#     sp = tuple(sp)  # convert to tuple
#     if sp == ():
#         return st  # nothing detected
#
#     if mode in [
#         "schmitt",
#         "Schmitt",
#     ]:  # normal operating mode is fixed voltage threshold, with hysterisis
#         for k in sp:
#             xx = xt[k - 1 : k + 1]
#             y = vma[k - 1 : k + 1]
#             if interpolate:
#                 m = (y[1] - y[0]) / dt  # local slope
#                 b = y[0] - (xx[0] * m)
#                 st = np.append(st, xx[1] + (thresh - b) / m)
#             else:
#                 if len(x) > 1:
#                     st = np.append(st, xx[1])
#
#     elif mode == "peak":
#         kpkw = int(peakwidth / dt)
#         z = (np.array(np.where(np.diff(spv) > 1)[0]) + 1).tolist()
#         #            print('z: ', z)
#         z.insert(0, 0)  # first element in spv is needed to get starting AP
#         for k in z:
#             zk = spv[k]
#             spk = np.argmax(vma[zk : zk + kpkw]) + zk  # find the peak position
#             xx = xt[spk - 1 : spk + 2]
#             y = vma[spk - 1 : spk + 2]
#             if interpolate:
#                 try:
#                     # mimic Igor FindPeak routine with B = 1
#                     m1 = (y[1] - y[0]) / dt  # local slope to left of peak
#                     b1 = y[0] - (xx[0] * m1)
#                     m2 = (y[2] - y[1]) / dt  # local slope to right of peak
#                     b2 = y[1] - (xx[1] * m2)
#                     mprime = (
#                         m2 - m1
#                     ) / dt  # find where slope goes to 0 by getting the line
#                     bprime = m2 - ((dt / 2.0) * mprime)
#                     st = np.append(st, -bprime / mprime + xx[1])
#                 except:
#                     continue
#             else:
#                 # print('utility: yere', x)
#                 if len(xx) > 1:
#                     st = np.append(st, xx[1])  # always save the first one
#
#     # clean spike times
#     # # st = clean_spiketimes(st, mindT=refract)
#     # print(("nspikes detected: ", len(st)), 'max spike time:', np.max(st))
#     # st2 = self.clean_spiketimes(st, mindT=refract)
#    # print(("nspikes detected after cleaning: ", len(st2)))


end
