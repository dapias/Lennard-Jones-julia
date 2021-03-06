# Geometric integrator for simulations in the canonical ensemble
## Lennard-Jones system

We provide the code that supports the numerical simulation reported in the manuscript [Geometric integrator for simulations in the canonical ensemble](https://aip.scitation.org/doi/abs/10.1063/1.4961506?journalCode=jcp). Note: the Logistic Thermostat was updated following the prescription stated in the paper [Ergodicity of One-dimensional Systems Coupled to the Logistic Thermostat](http://cmst.eu/articles/ergodicity-of-one-dimensional-systems-coupled-to-the-logistic-thermostat/). This change obeys to control the scale of the Logistic distribution instead of the mean in order to make the deterministic equations time-reversible.

It is organized as follows.

The ``src`` folder contains the files *MicrocanonicalSimulation.jl*, *ThermostatModel.jl*, *CanonicalSimulation.jl* and *LennardJonesModel.jl*. The first three files are generic and may be used to perform simulations in the microcanonical and canonical ensemble of any system with short-range interactions provided that the system is codified in a similar manner to that of Lennard Jones.

The module **CanonicalSimulation** is imported in the script *canonicalsimulation.jl* which takes certain values for the parameters (see below), performs the numerical integration and generates a *.hdf5* file which is saved in the ``HDF5`` folder with the name given by the user (asked by the script together with the kind of the thermostat).

The parameters reside in the file *parameters.yaml* and the users may modify it for their convenience.

Finally, the folder ``notebooks`` contains a *jupyter notebook* which illustrates how to manipulate the HDF5 data to generate the kind of figures displayed in the paper.

##Usage

Clone this repository with the following command in a UNIX terminal
```
~$ git clone https://github.com/dapias/Lennard-Jones-julia.git
```

Then move into the created folder and execute the script.  To do that you may proceed in one of the three different following ways

1. In a UNIX terminal type

 ```
 ~$ julia canonicalsimulation.jl
 ```
2. In a UNIX terminal execute julia as
 ```
 ~$ julia
 ```
This command opens the Julia [REPL](https://en.wikibooks.org/wiki/Introducing_Julia/The_REPL). And then type the following command
 ```
 julia> include("canonicalsimulation.jl")
 ```

3. Open a Jupyter Notebook and type in a cell (if the simulation is quite long, the browser might stop working!)
 ```
 include("canonicalsimulation.jl")
 ```

### Requirements

#### General
Julia. It may be downloaded from its [webpage](http://julialang.org/downloads/)

#### Particular
The following packages are needed for the adequate execution of the program

- HDF5
- YAML
- PyPlot (optional, it is used in the notebook)

To add a package type the following command in the Julia REPL.
```
julia> Pkg.add("PackageName")
```

### Miscellaneous
#### Gear Integration
The code for the Gear integration is in the branch *Gear*. Check the directory `src/` that contains the main code and the  directory `notebooks/` which contains a tutorial to run Gear integrations. 

To switch to the *Gear* branch: open a terminal, move to the directory *Lennard-Jones-julia* and then type

```
$ git checkout Gear
```

To go back to the master branch try
```
$ git checkout master
```
To just check the main code online click [here](https://github.com/dapias/Lennard-Jones-julia/blob/Gear/src/GearIntegration.jl)

### Author

**Diego Tapias** (Facultad de Ciencias, UNAM) diego.tapias@nucleares.unam.mx

*2016.*







