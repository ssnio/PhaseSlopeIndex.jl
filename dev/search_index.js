var documenterSearchIndex = {"docs":
[{"location":"examples/","page":"Examples","title":"Examples","text":"EditURL = \"https://github.com/ssnio/PhaseSlopeIndex.jl/blob/master/docs/literate/examples.jl\"","category":"page"},{"location":"examples/#Phase-Slope-Index-method","page":"Examples","title":"Phase Slope Index method","text":"","category":"section"},{"location":"examples/#Purpose","page":"Examples","title":"Purpose","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This is a walk-through notebook on Robustly Estimating the Flow Direction of Information in Complex Physical Systems paper, by Guido Nolte, Andreas Ziehe, Vadim V. Nikulin, Alois Schlögl, Nicole Krämer, Tom Brismar, and Klaus-Robert Müller, implemented in Julia-Language (please see http://doc.ml.tu-berlin.de/causality/ and Nolte et al. 2008).","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The notebook can be viewed here:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: binder)\n(Image: nbviewer)","category":"page"},{"location":"examples/#Acknowledgement","page":"Examples","title":"Acknowledgement","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This work was funded by the German Federal Ministry of Education and Research (BMBF) in the project ALICE III under grant ref. 01IS18049B.","category":"page"},{"location":"examples/#Load-packages","page":"Examples","title":"Load packages","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"using PhaseSlopeIndex\nusing Random: randperm\nusing MAT: matread\nusing Plots: plot, heatmap, cgrad\nusing DSP: blackman","category":"page"},{"location":"examples/#Data","page":"Examples","title":"Data","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"mixed_data contains four channels, where:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"channels 1 and 2 are i.i.d. uniform 0 1 noise\nchannel 3 is delayed (by 1 sample) channel 1\nchannel 4 is delayed (by 16 samples) channel 1 plus i.i.d. uniform 0 02 noise","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"# data generation\nT = Float64\nn_channels = 4  # number of channels\nn_samples = 2^16  # number of data points measured in each channel\nfs = 128  # sampling frequency\ntime_array = Array{T}(range(1; step=1 / fs, length=n_samples))\n\n# mixed data\nrange_c4 = 1:n_samples\nrange_c1 = range_c4 .+ 16\nrange_c3 = range_c1 .- 1\n\nrand_data = randn(T, (n_samples + 16, 1)) # uniform noise\ncause_source = rand_data[range_c1]  # channel 1\nrandom_source = randn(T, n_samples)  # channel 2, uniform noise\neffect_source = rand_data[range_c3]  #channel 3\nweak_effect = rand_data[range_c4] .- (randn(T, (n_samples, 1)) / 5) # channel 4\nmixed_data = hcat(cause_source, random_source, effect_source, weak_effect)\n\np1 = plot(\n    time_array[1:64],\n    mixed_data[1:64, :];\n    title=\"Mixed data\",\n    label=[\"Cause\" \"Random\" \"Effect\" \"Noisy Effect\"],\n    linecolor=[\"red\" \"green\" \"blue\" \"magenta\"],\n)\n\nplot(p1; layout=(1, 1), size=(800, 450))","category":"page"},{"location":"examples/#PSI","page":"Examples","title":"PSI","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"@doc data2psi","category":"page"},{"location":"examples/#Example-1","page":"Examples","title":"Example 1","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"PSI is calculated over all frequencies for segmented (seglen = 100) but continuous data (single epoch, nep = 1) and estimation of error using Bootstrap method for 256 resampling iterations (nboot=256). The default window function (Hanning window) is used.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"seglen = 100  # segment length\nnboot = 256  # number of bootstrap iterations\nmethod = \"bootstrap\"  # standard error estimation method\n\npsi, psi_se = data2psi(mixed_data, seglen; nboot=nboot, method=method)\n\np1 = heatmap(\n    psi;\n    ticks=false,\n    yflip=true,\n    yticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    xticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    color=cgrad(:bwr),\n    title=\"Phase Slope Index\",\n)\n\np2 = heatmap(\n    psi_se;\n    ticks=false,\n    yflip=true,\n    yticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    xticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    color=cgrad(:grays; rev=true),\n    title=\"PSI standard error\",\n)\n\nplot(p1, p2; layout=(1, 2), size=(720, 300))","category":"page"},{"location":"examples/#Example-2","page":"Examples","title":"Example 2","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"PSI is calculated over 3 frequency bands, for partitioned data to segments (seglen = 100) and epochs (eplen = 200), estimation of error using Jackknife method (default). The window function is set to blackman (imported from DSP.jl). The plots are for only one of the frequency ranges.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We normalize the PSI by dividing it by estimated standard error.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"seglen = 100  # segment length\neplen = 200  # epoch length\nmethod = \"jackknife\"  # standard error estimation method\n\n# three frequency bands\nfreqlist = [[5:1:10;] [6:1:11;] [7:1:12;]]\n\nsegave = true  # average across CS segments\nsubave = true  # subtract average across CS segments\ndetrend = true  # performs a 0th-order detrend across raw segments\nwindow = blackman  # blackman window function\n\npsi, psi_se = data2psi(\n    mixed_data,\n    seglen;\n    subave=subave,\n    segave=segave,\n    detrend=detrend,\n    freqlist=freqlist,\n    eplen=eplen,\n    method=method,\n    window=blackman,\n)\n\npsi_normed = psi ./ (psi_se .+ eps())\n\np1 = heatmap(\n    psi[:, :, 1];\n    ticks=false,\n    yflip=true,\n    yticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    xticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    color=cgrad(:bwr),\n    title=\"Phase Slope Index\",\n)\n\np2 = heatmap(\n    psi_se[:, :, 1];\n    ticks=false,\n    yflip=true,\n    yticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    xticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    color=cgrad(:grays; rev=true),\n    title=\"PSI standard error\",\n)\n\np3 = heatmap(\n    psi_normed[:, :, 1];\n    ticks=false,\n    yflip=true,\n    yticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    xticks=([1, 2, 3, 4], [\"Ch1\", \"Ch2\", \"Ch3\", \"Ch4\"]),\n    color=cgrad(:bwr),\n    title=\"Normalized PSI\",\n)\n\nplot(p1, p2, p3; layout=(1, 3), size=(1110, 270))","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#PhaseSlopeIndex.jl","page":"Home","title":"PhaseSlopeIndex.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is a JuliaLang implementation of Phase Slope Index method, proposed in \"Robustly Estimating the Flow Direction of Information in Complex Physical Systems\" , by Guido Nolte, Andreas Ziehe, Vadim V. Nikulin, Alois Schlögl, Nicole Krämer, Tom Brismar, and Klaus-Robert Müller (please see http://doc.ml.tu-berlin.de/causality/ and Nolte et al. 2008).","category":"page"},{"location":"#Acknowledgement","page":"Home","title":"Acknowledgement","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This work was funded by the German Federal Ministry of Education and Research (BMBF) in the project ALICE III under grant ref. 01IS18049B.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"index.md\", \"examples.md\"]","category":"page"},{"location":"#Functions","page":"Home","title":"Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The only exported function is data2psi:","category":"page"},{"location":"","page":"Home","title":"Home","text":"data2psi","category":"page"},{"location":"#PhaseSlopeIndex.data2psi","page":"Home","title":"PhaseSlopeIndex.data2psi","text":"data2psi(data, seglen [, segshift, eplen, freqlist, method,\n                         nboot, segave, subave, detrend])\n\ncalculates phase slope index (PSI)\n\nArguments\n\ndata::AbstractArray: NxM array for N data points in M channels\nseglen::Integer: segment length (determines the frequency resolution)\n\noptional arguments\n\nsegshift::Integer: number of bins by which neighboring segments are shifted (default=seglen/2)\neplen::Integer: length of epochs (if eplen=0, eplen is defaulted to number of samples)\nfreqlist::AbstractArray: 2D Array where each column is a frequency band (default is full range)\nmethod::String: standard deviation estimation method (default is \"jackknife\")\nnboot::Integer: number of bootstrap resamplings (default is 100)\nsegave::Bool: if true, average across CS segments (default is false)\nsubave::Bool: if true, subtract average across CS segments (default is false)\ndetrend::Bool: if true, performs a 0th-order detrend across raw segments (default is false)\nwindow::Function: window function with interval length as sole necessary argument (default is Hanning)\n\nReturns\n\npsi::AbstractArray: Phase Slope Index with shape (channel, channel, frequency bands)\npsi_se::AbstractArray: PSI estimated standard error with shape (channel, channel, frequency bands)\n\n\n\n\n\n","category":"function"}]
}