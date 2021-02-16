# PhaseSlopeIndex.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ssnio.github.io/PhaseSlopeIndex.jl/stable)
[![Master](https://img.shields.io/badge/docs-dev-blue.svg)](https://ssnio.github.io/PhaseSlopeIndex.jl/dev)
[![Build Status](https://github.com/ssnio/PhaseSlopeIndex.jl/workflows/CI/badge.svg)](https://ssnio.github.io/PhaseSlopeIndex.jl/actions)
[![Coverage](https://codecov.io/gh/ssnio/PhaseSlopeIndex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ssnio/PhaseSlopeIndex.jl)

This is a JuliaLang implementation of *Phase Slope Index* method, proposed in [*Robustly Estimating the Flow Direction of Information in Complex Physical Systems*](http://link.aps.org/abstract/PRL/v100/e234101). 
Please refer to http://doc.ml.tu-berlin.de/causality/ for more information.

### Usage
The only exported function is [`data2psi`](https://ssnio.github.io/PhaseSlopeIndex.jl/dev/#Functions).

For detailed examples, please refer to the [docs](https://ssnio.github.io/PhaseSlopeIndex.jl/dev/examples/).

### Citation
Please [cite the following paper](https://github.com/ssnio/PhaseSlopeIndex.jl/blob/master/citation.bib) if you use this code in published work:
> Nolte, G., Ziehe, A., Nikulin, V., Schlögl, A., Krämer, N., Brismar, T., & Müller, K.R. (2008), *[Robustly Estimating the Flow Direction of Information in Complex Physical Systems](http://link.aps.org/abstract/PRL/v100/e234101)*, Phys. Rev. Lett., 100, 234101. 

### Acknowledgement
This work was funded by the German Federal Ministry of Education and Research ([BMBF](https://www.bmbf.de/)) in the project ALICE III under grant ref. 01IS18049B.
