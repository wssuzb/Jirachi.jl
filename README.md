[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wssuzb.github.io/Jirachi.jl/dev/)
[![DOI](https://zenodo.org/badge/730485481.svg)](https://zenodo.org/doi/10.5281/zenodo.10428782)

# Jirachi - Julia project for analyzing variability in active galaxies harboring massive black holes

<img align="right" alt="jirachi" src="./test/fig/jirachi.jpeg" width="200" height="200"/>

<!-- - JIRACHI is a cute pokemon who always makes wishes come true !!! -->

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wssuzb.github.io/Jirachi.jl/stable/) -->

<!-- [![Build Status](https://github.com/wssuzb/Jirachi.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wssuzb/Jirachi.jl/actions/workflows/CI.yml?query=branch%3Amain) -->
<!-- [![Coverage](https://codecov.io/gh/wssuzb/Jirachi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/wssuzb/Jirachi.jl) -->

## About this project

Several functions can be found in this project to help us get better understand the central engine of AGNs:

- color variation
- structure function
- generate DRW

A `PYTHON` version for calculating color variation can be found at [pyColorVariation.py](https://github.com/wssuzb/Jirachi.jl/blob/main/py/pyColorVariation.py).

Read the [full documentation here.](https://wssuzb.github.io/Jirachi.jl/dev/)

## Install

```julia
julia> # Press the key "]"

(@v1.9) pkg> add https://github.com/wssuzb/Jirachi.jl.git
julia> using Jirachi
```

or alternatively, you can also download this `Jirachi.jl-main.zip` file, and load it by
```julia
julia> push!(LOAD_PATH, path + "Jirachi.jl-main")
julia> using Jirachi
```
you just need to modify the `path`, where the `Jirachi` is downloaded, and enjoy youself!

<!-- ## -->

<!-- - structure function


<img align="center" alt="sf" src="./test/fig/plot_sf.svg" width="500" height="500"/> -->


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
