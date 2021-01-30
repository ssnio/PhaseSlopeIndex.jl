isdefined(Base, :__precompile__) && __precompile__()

module PhaseSlopeIndex

using Statistics: mean, std
using FFTW: fft!
using DSP: hanning
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
dropmean(X, d) = dropdims(mean(X, dims=d), dims=d)


"""
    squeeze(X::AbstractArray)

removing singleton dimensions
"""
function squeeze(X::AbstractArray)
    keepd = ()
    for (i, d) in enumerate(size(X))
        if d != 1; keepd = (keepd..., d) end
    end
    reshape(X, keepd)
end


"""
    data2para(data, seglen, segshift, eplen, freqlist,
              method, nboot, segave, subave, detrend)

Extracts and builds a named tuple of parameters.

# Arguments
- `data::AbstractArray`: NxM array for N data points in M channels.
- `seglen::Integer`: segment length (determinds the frequency resolution).
- `segshift::Integer`: number of bins by which neighboring segments are shifted.
  e.g. segshift=seglen/2 makes overlapping segments
- `eplen::Integer`: length of epochs
- `freqlist::AbstractArray`: 2D Array where each column is a frequency band
- `method::String`: standard deviation estimation method
- `nboot::Integer`: number of bootstrap resamplings
- `segave::Bool`: if true, average across segments for CS calculation
- `subave::Bool`: if true, subtract average across segments and epochs for CS calculation
  (**IMPORTANT**: For just one epoch (e.g. for continuous data) set subave = false)
- `detrend::Bool`: if true, performes a linear detrend across segments

# Returns
- `parameters::NamedTuple`: a named tuple of parameters

"""
function data2para(data::AbstractArray,
                   seglen::Integer,
                   segshift::Integer,
                   eplen::Integer,
                   freqlist::AbstractArray,
                   method::String,
                   nboot::Integer,
                   segave::Bool,
                   subave::Bool,
                   detrend::Bool)

    # We would like to avoid transpose and copying the data!
    if size(data, 1) < size(data, 2)
        @info "data is transposed to (#samples, #channels) "
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
    nep = int(nsamples/eplen)  # number of epochs

    if segshift == 0 ; segshift = int(seglen/2) end
    nseg = int((eplen - seglen) / segshift) + 1

    # size(freqlist) = (freqs, nfbands)
    if length(freqlist) == 0; freqlist = reshape(Array(1:int(seglen/2)), (:, 1)) end
    maxfreq = maximum(freqlist)  # TODO: maximum(freqlist, dims=1)
    nfbands = size(freqlist, 2)  # multiple frequency bands

    # we use named tuples to book the parameters
    parameters = (data=data,
                  nsamples=nsamples,
                  nchan=nchannels,
                  eplen=eplen,
                  nep=nep,
                  method=method,
                  nboot=nboot,
                  seglen=seglen,
                  segshift=segshift,
                  nseg=nseg,
                  freqlist=freqlist,
                  maxfreq=maxfreq,
                  nfbands=nfbands,
                  segave=segave,
                  subave=subave,
                  detrend=detrend)

    return parameters
end


"""
    make_eposeg(data, seglen, nep, nseg, nchan, segshift)

Partitioning data into epochs and segments

# Arguments
- `data::AbstractArray`: NxM array for N data points in M channels.
- `seglen::Integer`: segment length
- `nep::Integer`: number of epochs
- `nseg::Integer`: number of segments
- `nchan::Integer`: number of channels
- `segshift::Integer`: number of bins by which neighboring segments are shifted.

# Returns
- `epseg::AbstractArray`: partitioned data into shape (seglen, nep, nseg, nchan)

**Note**: The returned `epseg` may have more data entries than the input data.

"""
function make_eposeg(data::AbstractArray,
                     seglen::Integer,
                     eplen::Integer,
                     nep::Integer,
                     nseg::Integer,
                     nchan::Integer,
                     segshift::Integer)::AbstractArray

    # preallocation
    epseg = Array{Complex{Float64}}(undef, seglen, nep, nseg, nchan)

    for (i, e) in zip(1:nep, 1:eplen:nep*eplen)
        for (j, s) in zip(1:nseg, 1:segshift:nseg*seglen)
            @inbounds @views epseg[:, i, j, :] = data[e:e+eplen-1, :][s:s+seglen-1, :]
        end
    end

    return epseg
end


"""
    detrend!(data)
(in place) Linear detrend of signals along first axis

# Arguments
- `data::AbstractArray`: N-dim array where signal is stored in column-major order

**Note**: shape of data must be (signal length, ...)

"""
function detrend!(data::AbstractArray)

    original_shape = size(data)
    nsamp = size(data, 1)  # number of samples

    A = Array{Float64}([ones(nsamp) Array(1:nsamp)])

    data = reshape(data, (nsamp, :))  # reshaping data
    data .-= A * (A\data)
    reshape(data, original_shape)
end


"""
    cs2ps(cs)
Cross Spectra to Phase Slope

# Arguments
- `cs::AbstractArray`: the cross spectral array with size (seglen, :, nchan, nchan)

# Returns
- phase slope index (AbstractArray) as Eq. 3 of the reference paper

**Note**: the frequency resolution is assumed to be the resolution of freq. band!

"""
function cs2ps(cs::AbstractArray)

    # size(cs) =!= (seglen, nchan, nchan)
    # we don't use df but we slice the fband

    # complex coherency
    @einsum cx_coh[f, i, j] := cs[f, i, j] / sqrt(cs[f, i, i] * cs[f, j, j])

    # phase slope (Eq. 3)
    @views imag.(sum(conj(cx_coh[1:end-1, :, :]) .* cx_coh[2:end, :, :], dims=1))
end


"""
    data2ps(data)
Epoched data to Cross Spectra

# Arguments
- `data::AbstractArray`: segmented and epoched data in shape of (maxfreq, nep, nseg, nchan)

# Return
- `cs::AbstractArray{Complex}`: cross spectral as eq. 2 of reference paper.

"""
function data2cs(data::AbstractArray)
    # cs: cross-spectral matrix

    # Eq. 2 size(csepseg) = (maxfreq, nep, nseg, nchan, nchan)
    @einsum cs[i, j, k, m, n] := data[i, j, k, m] * conj(data[i, j, k, n])
end


"""
    data2cs(data, nboot, method)
Epoched data to Cross Spectra

# Arguments
- `data::AbstractArray`: cross spectral data in shape of (maxfreq, nep, nseg, nchan, nchan)
- 'nboot::Integer`: number of bootstrap epochs
- `method::String`: resampling method ("jackknife" or "bootstrap")

# Return
- `psi_rs::AbstractArray`: resampling PSI in shape of (nchan, nchan, nboot or nep)

"""
function cs2ps_std(data::AbstractArray,
                   nboot::Integer,
                   method::String)

    maxfreq, nep, nseg, nchan, nchan = size(data)

    if method == "jackknife"
        psi_rs = Array{Float64}(undef, nchan, nchan, nep)
        for z in 1:nep
            sliced_cs = view(data, :, (1:nep .!=z), :, :, :)
            @fastmath psi_rs[:, :, z] = cs2ps(dropmean(sliced_cs, (2, 3)))
        end
    elseif method == "bootstrap"
        psi_rs = Array{Float64}(undef, nchan, nchan, nboot)
        for z in 1:nboot
            sliced_cs = view(data, :, rand(1:nep, nep), :, :, :)
            @fastmath psi_rs[:, :, z] = cs2ps(dropmean(sliced_cs, (2, 3)))
        end
    end
    return psi_rs
end



"""
    data2psi(data, seglen [, segshift, eplen, freqlist, method,
                             nboot, segave, subave, detrend])
calculates phase slope index (PSI)

# Arguments
- `data::AbstractArray`: NxM array for N data points in M channels.
- `seglen::Integer`: segment length (determinds the frequency resolution).
- `segshift::Integer`: number of bins by which neighboring segments are shifted.
e.g. segshift=seglen/2 makes overlapping segments
- `eplen::Integer`: length of epochs
- `freqlist::AbstractArray`: 2D Array where each column is a frequency band
- `method::String`: standard deviation estimation method
- `nboot::Integer`: number of bootstrap resamplings
- `segave::Bool`: if true, average across segments for CS calculation
- `subave::Bool`: if true, subtract average across segments and epochs for CS calculation
(**IMPORTANT**: For just one epoch (e.g. for continuous data) set subave = false)
- `detrend::Bool`: if true, performes a linear detrend across segments

# Returns
- `psi::AbstractArray`: channel x channel PSI
- `psi_se::AbstractArray`: channel x channel PSI standard error
- `psi_normed::AbstractArray`: normalized PSI by sqrt(nrsmpl)*psi_se
"""
function data2psi(data::AbstractArray,
                 seglen::Integer;
                 segshift::Integer=0,
                 eplen::Integer=0,
                 freqlist::AbstractArray=Int64[],
                 method::String="bootstrap",
                 nboot::Integer=100,
                 segave::Bool=false,
                 subave::Bool=false,
                 detrend::Bool=false)

    para = data2para(data, seglen, segshift, eplen, freqlist,
                     method, nboot, segave, subave, detrend)

    eposeg_data = make_eposeg(para.data, para.seglen, para.eplen, para.nep,
                              para.nseg, para.nchan, para.segshift)

    if para.detrend; detrend!(eposeg_data) end

    eposeg_data .*= hanning(para.seglen)

    fft!(eposeg_data, 1)

    if para.method == "jackknife"
        nrsmpl = para.nep
    elseif para.method == "bootstrap"
        nrsmpl = para.nboot
    end

    # preallocation
    psi = Array{Float64}(undef, para.nfbands, para.nchan, para.nchan)
    psi_rs = Array{Float64}(undef, para.nfbands, para.nchan, para.nchan, nrsmpl)

    for (f, fband) in enumerate(eachrow(para.freqlist'))
        cs_f = data2cs(view(eposeg_data, fband, :, :, :))
        psi[f, :, :] = cs2ps(dropmean(cs_f, (2, 3)))
        psi_rs[f, :, :, :] = cs2ps_std(cs_f, para.nboot, para.method)
    end

    psi_se = squeeze(std(psi_rs, mean=psi, dims=4))
    psi = squeeze(psi)
    psi_normed = psi ./ (sqrt(para.nep) .* psi_se)

    return psi, psi_se, psi_normed
end


end  # module
