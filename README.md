[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wssuzb.github.io/Jirachi.jl/dev/)
[![DOI](https://zenodo.org/badge/730485481.svg)](https://zenodo.org/doi/10.5281/zenodo.10428782)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wssuzb.github.io/Jirachi.jl/stable/)
[![Build Status](https://github.com/wssuzb/Jirachi.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wssuzb/Jirachi.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/wssuzb/Jirachi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/wssuzb/Jirachi.jl)

# Jirachi - Julia project for analyzing variability in active galaxies harboring massive black holes

<img align="right" alt="jirachi" src="./test/fig/jirachi.jpeg" width="200" height="200"/>

## About this project

Several functions can be found in this project to help us get a better understand of the central engine of AGNs:

- timescale-dependent color variation calculation.
- structure function calculation and fits.
- generate stochastic process light curves.
- fractional variability.
- ...

A `PYTHON` version for calculating color variation can be found at [run_demo.py](https://github.com/wssuzb/Jirachi.jl/blob/main/py/run_demo.py).

<img align="center" alt="jirachi" src="./py/mag_flux_new_sun14_i_z.png" width="700" height="600"/>


Read the [full documentation here.](https://wssuzb.github.io/Jirachi.jl/dev/)

TODO:

- [ ] reprocessing model
- [ ] structure function calculation (see kozlowski+16)

## Install

```julia
julia> # Press the key "]"

(@v1.9) pkg> add https://github.com/wssuzb/Jirachi.jl.git
julia> using Jirachi
```

or alternatively, you can also download this `Jirachi.jl-main.zip` file, and load it by
```julia
julia> push!(LOAD_PATH, "'~/where/you/download/the/package/Jirachi")
julia> using Jirachi
```
enjoy!


## Citing

If this project makes you life easier, pls. cites this code as below:

```bib
@software{jirachi,
  author       = {Su, Zhen-Bo},
  title        = {{Julia project for analyzing variability in active 
                   galaxies harboring massive black holes}},
  month        = dec,
  year         = 2023,
  publisher    = {Zenodo},
  version      = {v0.0.1},
  doi          = {10.5281/zenodo.10428783},
  url          = {https://doi.org/10.5281/zenodo.10428783}
}
```

also see [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
