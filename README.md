#Taylor series based integrator for contact hamiltonian flows

We provide the code that supports the numerical simulation reported in the manuscript [A thermostat algorithm generating target ensembles](http://arxiv.org/abs/1510.03942).

The code is based on the package  [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl).

It is organized as follows.

The ``src`` folder contains the file  *ContactIntegrator.jl*. This file defines the module **ContactIntegrator** that exports the main function **contacthointegration!** that performs the core of the simulation.

The module is imported in the script *contact_taylor.jl* which takes certain values for the parameters (see below), performs the numerical integration and generates a *.hdf5* file which is saved in the ``HDF5`` folder with the name given by the user (asked by the script).

The parameters reside in the file *parameters.yaml* and the users may modify it for their convenience.

Finally, the folder ``notebooks`` contains a *jupyter notebook* which illustrates how to manipulate the HDF5 data to generate the kind of figures displayed in the paper.

##Usage

Clone this repository with the following command in a UNIX terminal
```
~$ git clone https://github.com/dapias/ContactFlowsTaylor.git
```

Then move into the created folder and execute the script.  To do that you may proceed in one of the three different following ways

1. In a UNIX terminal type

 ```
 ~$ julia contact_taylor.jl
 ```
2. In a UNIX terminal execute julia as
 ```
 ~$ julia
 ```
This command opens the Julia [REPL](https://en.wikibooks.org/wiki/Introducing_Julia/The_REPL). And then type the following command
 ```
 julia> include("contact_taylor.jl")
 ```

3. Open a Jupyter Notebook and type in a cell
 ```
 include("contact_taylor.jl")
 ```

### Requirements

#### General
Julia. It may be downloaded from its [webpage](http://julialang.org/downloads/)

#### Particular
The following packages are needed for the adequate execution of the program

- TaylorSeries
- HDF5
- YAML
- PyPlot (optional, it is used in the notebook)

To add a package type the following command in the Julia REPL.
```
julia> Pkg.add("PackageName")
```
###Miscellaneous (Runge-Kutta)

In the folder RungeKutta, the same field was integrated using the  4th order adaptive Runge-Kutta adaptive solver with the [Dormand-Price](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method) method implemented in the Package [ODE.jl](https://github.com/JuliaLang/ODE.jl). It can be checked that the claim of the manuscript does not depend on the method of integration by analyzing the results of the Runge-Kutta integration. However, for this procedure there is a systematic numeric error that even though it is small it tends to grow up with time. It can be analyzed following the evolution of the invariant quantity.

###Authors

**Diego Tapias** (Facultad de Ciencias, UNAM) diego.tapias@nucleares.unam.mx

**Alessandro Bravetti** (Instituto de Ciencias Nucleares, UNAM) bravetti@correo.nucleares.unam.mx

*2015.*







