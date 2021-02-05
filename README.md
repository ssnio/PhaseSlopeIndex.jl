# PhaseSlopeIndex.jl

This is a JuliaLang implementation of *Phase Slope Index* method, proposed in "*Robustly Estimating the Flow Direction of Information in Complex Physical Systems*" , by *Guido Nolte, Andreas Ziehe, Vadim V. Nikulin, Alois Schlögl, Nicole Krämer, Tom Brismar, and Klaus-Robert Müller* (please see http://doc.ml.tu-berlin.de/causality/ and [Nolte et al. 2008](http://link.aps.org/abstract/PRL/v100/e234101)).

### Acknowledgement
This work was funded by the German Federal Ministry of Education and Research [(BMBF)](https://www.bmbf.de/) in the project ALICE III under grant ref. 01IS18049B.

### Usage

The only exported function is `data2psi`:

```julia
psi, psi_se = function data2psi(
      data::AbstractArray,
      seglen::Integer;
      segshift::Integer = 0,
      eplen::Integer = 0,
      freqlist::AbstractArray = Int64[],
      method::String = "bootstrap",
      nboot::Integer = 100,
      segave::Bool = false,
      subave::Bool = false,
      detrend::Bool = false,
      window)
```

### Arguments
- `data::AbstractArray`: NxM array for N data points in M channels.
- `seglen::Integer`: segment length (determinds the frequency resolution).
- `segshift::Integer`: number of bins by which neighboring segments are shifted.
 e.g. segshift=seglen/2 makes overlapping segments
- `eplen::Integer`: length of epochs
- `freqlist::AbstractArray`: 2D Array where each column is a frequency band
- `method::String`: standard deviation estimation method
- `nboot::Integer`: number of bootstrap resamplings
- `segave::Bool`: if true, average across segments for CS calculation
- `subave::Bool`: if true, subtract average across segments for CS calculation
(For single epoch (e.g. for continuous data) set subave = false)
- `detrend::Bool`: if true, performes a linear detrend across segments
- `window`: the window function with interval length as sole necessary argument

### Returns
- `psi::AbstractArray`: channel x channel PSI
- `psi_se::AbstractArray`: channel x channel PSI estimated standard error
