# PhaseSlopeIndex.jl

This is a JuliaLang implementation of *Phase Slope Index* method, proposed in "*Robustly Estimating the Flow Direction of Information in Complex Physical Systems*" , by *Guido Nolte, Andreas Ziehe, Vadim V. Nikulin, Alois Schlögl, Nicole Krämer, Tom Brismar, and Klaus-Robert Müller* (please see http://doc.ml.tu-berlin.de/causality/ and [Nolte et al. 2008](http://link.aps.org/abstract/PRL/v100/e234101)).

### Acknowledgement
This work was funded by the German Federal Ministry of Education and Research [(BMBF)](https://www.bmbf.de/) in the project ALICE III under grant ref. 01IS18049B.

### Usage

For detailed examples, please see the [PhaseSlopeIndex.jl.ipynb](https://github.com/ssnio/PhaseSlopeIndex.jl/blob/master/notebooks/PhaseSlopeIndex.jl.ipynb).

The only exported function is `data2psi`:

```julia
psi, psi_se = function data2psi(
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
```

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
