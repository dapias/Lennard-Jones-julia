push!(LOAD_PATH,"src/")

import YAML
using CanonicalSimulation
using HDF5

println("Type the kind of the thermostat (Gaussian, Logistic or Quartic) \n Case sensitive")
thermostat = string(readline(STDIN))
thermoname = thermostat[1:end-1]

thermolist = ["Gaussian", "Logistic", "Quartic"]

while !(thermoname in thermolist)
  println("The thermostat kind you typed is not in our database. Try one of the following: \n Gaussian, Logistic or Quartic or check the spelling")
  thermostat = string(readline(STDIN))
  thermoname = thermostat[1:end-1]
end



println("Type the output filename (without a format):")
name = string(readline(STDIN))
filename = name[1:end-1]

while isfile("./data/$(thermoname)/HDF5/$filename.hdf5")
  println("The filename typed already exists in the HDF5 folder. Try another one:")
  filename = string(readline(STDIN))
end

parameters = YAML.load(open("parameterscanonical.yaml"))

N = parameters["N"]
rho = parameters["rho"]
T = parameters["T"]
dt = parameters["dt"]
runtime = parameters["runtime"]
Q = parameters["Q"]

file = h5open("./data/$(thermoname)/HDF5/$filename.hdf5", "w")

attrs(file)["Thermostat"] = thermoname
attrs(file)["N"] = N
attrs(file)["T"] = T
attrs(file)["deltat"] = dt
attrs(file)["Q"] = Q

results =  CanonicalSimulation.run(runtime, rho, dt, T, N, Q, thermoname)


time = results[1]
energy = results[2]
kinetic = results[3]
potential = results[4]
temperature = results[5]
invariant = results[6]
atoms = results[7]   ###For the moment, I don't save it in HDF5 file
peta = results[8]
etas = results[9]
vrandatom = results[10]

file["t"] = time
attrs(file)["runtime"] = time[end]
file["E"] = energy
file["K"] = kinetic
file["U"] = potential
file["T"] = temperature
file["invariant"] = invariant
file["peta"] = peta
file["eta"] = etas
file["vrandatom"] = vrandatom


close(file)

println("The simulation was succesfully done. \nOutput in .
/data/$(thermoname)/HDF5/$filename.hdf5")
