# CardioModels.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg) [![codecov.io](http://codecov.io/github/mapi1/CardioModels.jl/coverage.svg?branch=master)](http://codecov.io/github/mapi1/CardioModels.jl?branch=master) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mapi1.github.io/CardioModels.jl/dev)

This little package gives some models of cardiovascular variability, for experimenting and learning.

# Installation

To install this package type the following into the REPL
```julia
using Pkg
Pkg.add("https://github.com/mapi1/CardioModels.jl.git")
```

# Examples

Example notebooks using [Pluto.jl](https://github.com/fonsp/Pluto.jl) are provided in the folder `/examples` showing basic usage examples. The folder contains its own `Project.toml`. By running Julia in the example folder and using

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

all needed packages will be installed.


# Models
Below is a minimal introduction of the models, for more details refer to my thesis, the docs or the original publications.
## DeBoer-model

## Karemaker-model

## Baselli-model
