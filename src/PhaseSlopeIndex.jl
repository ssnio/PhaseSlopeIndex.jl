
module PhaseSlopeIndex

using Statistics: mean, std
using FFTW: fft
using Einsum

# Exports
#---
export data2psi


"""
    int(x) = floor(Int64, x)
"""
int(x) = floor(Int64, x)


"""
    dropmean(X, d) = dropdims(mean(X, dims=d), dims=d)
"""
dropmean(X, d) = dropdims(mean(X, dims = d), dims = d)


"""
    squeeze(X::AbstractArray)

removing singleton dimensions
"""
function squeeze(X::AbstractArray)
    keepd = ()
    for (i, d) in enumerate(size(X))
        if d != 1 ; keepd = (keepd..., d) end
    end
    reshape(X, keepd)
end


"""
    detrend!(data, n)
(in place) Linear detrend of signals along first axis

### Arguments
- `data::AbstractArray`: N-dim array where signal is in column-major order
- `n::Integer`: n=0 subtracts mean from data, n=1 removes linear trend

**Note**: shape of data must be (signal length, ...)

"""
function detrend!(data::AbstractArray, n::Integer)

    original_shape = size(data)
    nsamp = size(data, 1)  # number of samples

    A = Array{Float64}([ones(nsamp) Array(1:nsamp)])

    data = reshape(data, (nsamp, :))  # reshaping data
    if n == 0
        data .-= mean(data, dims = 1)
    elseif n == 1
        data .-= A * (A \ data)
    end
    reshape(data, original_shape)
end


"""
    window = hanning_fun(N)

Hanning window similar to MATLAB implementation
"""
function hanning_fun(N::Integer)
    x = [range(0.0, 1.0, length = N + 2);]
    window = map(x -> 0.5 * (1 - cospi(2 * x)), x)
    return window[2:end-1]
end


"""
    data2para(data, seglen, segshift, eplen, freqlist, method, nboot)

Extracts and builds a named tuple of parameters.

### Arguments
- `data::AbstractArray`: NxM array for N data points in M channels.
- `seglen::Integer`: segment length (determines the frequency resolution).
- `segshift::Integer`: number of bins by which neighboring segments are shifted.
e.g. segshift=seglen/2 makes overlapping segments
- `eplen::Integer`: length of epochs
- `freqlist::AbstractArray`: 2D Array where each column is a frequency band
- `method::String`: standard deviation (error) estimation method
- `subave::Bool`: if true, subtract average from CS segments (for continuous data, subave = false)
- `nboot::Integer`: number of bootstrap resamplings

### Returns
- `parameters::NamedTuple`: a named tuple of parameters

"""
function data2para(
    data::AbstractArray,
    seglen::Integer,
    segshift::Integer,
    eplen::Integer,
    freqlist::AbstractArray,
    method::String,
    subave::Bool,
    nboot::Integer)

    # We would like to avoid transpose and copying the data!
    if size(data, 1) < size(data, 2)
        @info "data is transposed to (#samples, #channels)"
        data = reshape(data, size(data, 2), size(data, 1))
    end

    # number of samples per channel and number of channels
    nsamples, nchannels = size(data)

    # jackknife check
    if eplen == 0 && method == "jackknife"
        @warn "Epoch length = 0 => Bootstrap method will be used."
        method = "bootstrap"
    end

    # if eplen = nsamples: continuous recording
    if eplen == 0 ; eplen = nsamples end
    nep = int(nsamples / eplen)  # number of epochs

    if nep == 1 && subave == true
        subave = false
        @warn "subave is set to false for continuous data"
    end

    if segshift == 0 ; segshift = int(seglen / 2) end
    nseg = int((eplen - seglen) / segshift) + 1

    # size(freqlist) = (freqs, nfbands)
    if length(freqlist) == 0
        freqlist = reshape(Array(1:int(seglen / 2)+1), (:, 1))
    elseif ndims(freqlist) == 1
        freqlist = reshape(freqlist, :, 1)
    end
    if size(freqlist, 1) < size(freqlist, 2)
        @info "freqlist is transposed to (#freq, #nfbands)"
        freqlist = freqlist'
    end
    maxfreq = maximum(freqlist)  # TODO: maximum(freqlist, dims=1)
    nfbands = size(freqlist, 2)  # multiple frequency bands

    if method == "jackknife"
        n_resample = nep
    elseif method == "bootstrap"
        n_resample = nboot
    else
        throw(method * " is not a supported error estimation method!" *
                       " Please choose from [jackknife, bootstrap] methods.")
    end

    # we use named tuples to book the parameters
    parameters = (
        data = data,
        nsamples = nsamples,
        nchan = nchannels,
        eplen = eplen,
        nep = nep,
        method = method,
        subave = subave,
        segshift = segshift,
        nseg = nseg,
        freqlist = freqlist,
        maxfreq = maxfreq,
        nfbands = nfbands,
        n_resample = n_resample,
    )

    return parameters
end


"""
    make_eposeg(data, seglen, nep, nseg, nchan, segshift)

Partitioning data into epochs and segments

### Arguments
- `data::AbstractArray`: NxM array for N data points in M channels.
- `seglen::Integer`: segment length
- `nep::Integer`: number of epochs
- `nseg::Integer`: number of segments
- `nchan::Integer`: number of channels
- `segshift::Integer`: number of bins by which neighboring segments are shifted.

### Returns
- `epseg::AbstractArray`: partitioned data into shape (seglen, nep, nseg, nchan)

**Note**: returned `epseg` may have more data entries than input data.

"""
function make_eposeg(
    data::AbstractArray,
    seglen::Integer,
    eplen::Integer,
    nep::Integer,
    nseg::Integer,
    nchan::Integer,
    segshift::Integer)::AbstractArray

    # preallocation
    epseg = Array{Float64}(undef, seglen, nep, nseg, nchan)

    for (i, e) in zip(1:nep, 1:eplen:nep*eplen)
        for (j, s) in zip(1:nseg, 1:segshift:nseg*seglen)
            @views epseg[:, i, j, :] = data[e:e+eplen-1, :][s:s+seglen-1, :]
        end
    end

    return epseg
end


"""
    cs2ps(cs)
Cross Spectra to Phase Slope

### Arguments
- `cs::AbstractArray`: Cross Spectral array with size (seglen, :, nchan, nchan)

### Returns
- phase slope index (AbstractArray) as Eq. 3 of the reference paper

**Note**: frequency resolution is assumed to be the resolution of freq. band!

"""
function cs2ps(cs::AbstractArray)

    # size(cs) = (seglen, nchan, nchan)
    # we don't use df but we slice the fband

    # complex coherency
    @einsum coh[f, i, j] := cs[f, i, j] / sqrt(cs[f, i, i] * cs[f, j, j])

    # phase slope (Eq. 3)
    @views imag.(sum(conj(coh[1:end-1, :, :]) .* coh[2:end, :, :], dims = 1))
end


"""
    data2ps(data)
Epoched segmented data to Cross Spectra

### Arguments
- `data::AbstractArray`: Segmented data of shape (maxfreq, nep, nseg, nchan)

### Return
- `cs::AbstractArray{Complex}`: Cross Spectral as eq. 2 of reference paper.

"""
function data2cs(data::AbstractArray)
    # cs: cross-spectral matrix

    # Eq. 2 size(csepseg) = (maxfreq, nep, nseg, nchan, nchan)
    @einsum cs[f, e, s, i, j] := data[f, e, s, i] * conj(data[f, e, s, j])
end


"""
    _cs_ = cs2cs_(data, cs, fband, nep, segave, subave, method)

preparing Cross Spectra for Phase Slope by segment averaging and subtraction

### Arguments
- `data::AbstractArray`: Fourier-transformed detrended epoched segmented data.
- `cs::AbstractArray`: Cross Spectra of data
- `fband::AbstractArray`: 1D array of frequency range for PSI calculation
- `nep::Integer`: number of epochs
- `segave::Bool`: if true, averages across CS segments
- `subave::Bool`: if true, subtract average across CS segments
- `method::String`: standard deviation estimation method

### Returns
- `_cs_::AbstractArray`: segment averaged and subtracted Cross Spectra

"""
function cs2cs_(
    data::AbstractArray,
    cs::AbstractArray,
    fband::AbstractArray,
    nep::Integer,
    segave::Bool,
    subave::Bool,
    method::String)

    if segave
        if method == "bootstrap"
            randboot = rand(1:nep, nep)
            cs_ = dropmean(view(cs, :, randboot, :, :, :), (2, 3))
            av_ = dropmean(view(data, fband, randboot, :, :), (2, 3))
        elseif method == "psi"
            cs_ = dropmean(cs, (2, 3))
            av_ = dropmean(view(data, fband, :, :, :), (2, 3))
        elseif method == "jackknife"
            cs_ = dropmean(cs, 3)
            av_ = dropmean(view(data, fband, :, :, :), 3)
        end
    else
        if method == "bootstrap"
            randboot = rand(1:nep, nep)
            cs_ = dropmean(view(cs, :, randboot, 1, :, :), 2)
            av_ = dropmean(view(data, fband, randboot, 1, :), 2)
        elseif method == "psi"
            cs_ = dropmean(view(cs, :, :, 1, :, :), 2)
            av_ = dropmean(view(data, fband, :, 1, :), 2)
        elseif method == "jackknife"
            cs_ = view(cs, :, :, 1, :, :)
            av_ = view(data, fband, :, 1, :)
        end
    end

    if subave
        if method == "bootstrap" || method == "psi"
            @einsum av_ij[f, i, j] := av_[f, i] * conj(av_[f, j])
        elseif method == "jackknife"
            @einsum av_ij[f, e, i, j] := av_[f, e, i] * conj(av_[f, e, j])
        end
        _cs_ = squeeze(cs_ - av_ij)
    else
        if method == "bootstrap" || method == "psi" || method == "jackknife"
            _cs_ = squeeze(cs_)
        end
    end

    return _cs_
end


"""
    data2psi(data, seglen [, segshift, eplen, freqlist, method,
                             nboot, segave, subave, detrend])
calculates phase slope index (PSI)

### Arguments
- `data::AbstractArray`: NxM array for N data points in M channels
- `seglen::Integer`: segment length (determines the frequency resolution)

*optional arguments*
- `segshift::Integer`: number of bins by which neighboring segments are shifted (default=seglen/2)
- `eplen::Integer`: length of epochs (if eplen=0, eplen is defaulted to number of samples)
- `freqlist::AbstractArray`: 2D Array where each column is a frequency band (default is full range)
- `method::String`: standard deviation estimation method (default is "jackknife")
- `nboot::Integer`: number of bootstrap resamplings (default is 100)
- `segave::Bool`: if true, average across CS segments (default is false)
- `subave::Bool`: if true, subtract average across CS segments (default is false)
- `detrend::Bool`: if true, performs a 0th-order detrend across raw segments (default is false)
- `window::Function`: window function with interval length as sole necessary argument (default is Hanning)

### Returns
- `psi::AbstractArray`: Phase Slope Index with shape (channel, channel, frequency bands)
- `psi_se::AbstractArray`: PSI estimated standard error with shape (channel, channel, frequency bands)

"""
function data2psi(
    data::AbstractArray,
    seglen::Integer;
    segshift::Integer = 0,
    eplen::Integer = 0,
    freqlist::AbstractArray = Int64[],
    method::String = "jackknife",
    nboot::Integer = 100,
    segave::Bool = false,
    subave::Bool = false,
    detrend::Bool = false,
    window::Function = hanning_fun)

    (
    data, nsamples, nchan, eplen, nep, method, subave,
    segshift, nseg, freqlist, maxfreq, nfbands, n_resample,
    ) = data2para(data, seglen, segshift, eplen, freqlist, method, subave, nboot)

    eposeg = make_eposeg(data, seglen, eplen, nep, nseg, nchan, segshift)

    if detrend
        detrend!(eposeg, 0)
    end

    if isa(window, Function)
        window = window(seglen)
    else
        throw(window * " is not a Function!")
    end
    eposeg .*= window

    eposeg = view(fft(eposeg, 1), 1:maxfreq, :, :, :)

    # preallocation
    psi = Array{Float64}(undef, nchan, nchan, nfbands)
    psi_est = Array{Float64}(undef, nchan, nchan, nfbands, n_resample)
    for (f, fband) in enumerate(eachrow(freqlist'))
        cs_full = data2cs(view(eposeg, fband, :, :, :))

        cs_psi = cs2cs_(eposeg, cs_full, fband, nep, segave, subave, "psi")
        psi[:, :, f] = cs2ps(cs_psi)

        if method == "jackknife"
            cs_jack = cs2cs_(eposeg, cs_full, fband, nep, segave, subave, "jackknife")
            for e = 1:nep
                cs_jack_se = (nep * cs_psi - view(cs_jack, :, e, :, :)) / (nep + 1)
                psi_est[:, :, f, e] = cs2ps(cs_jack_se)
            end

        elseif method == "bootstrap"
            for e = 1:nboot
                cs_boot = cs2cs_(eposeg, cs_full, fband, nep, segave, subave, "bootstrap")
                psi_est[:, :, f, e] = cs2ps(cs_boot)
            end
        end
    end

    psi = squeeze(psi)
    if method == "jackknife"
        psi_se = sqrt(nep) * squeeze(std(psi_est, corrected = true, dims = 4))
    elseif method == "bootstrap"
        psi_se = squeeze(std(psi_est, corrected = true, dims = 4))
    end

    return psi, psi_se
end


end  # module
