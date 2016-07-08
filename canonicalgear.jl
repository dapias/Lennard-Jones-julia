push!(LOAD_PATH,"src/")

import YAML
using GearIntegration
using HDF5

try
  mkdir("./data/")
  mkdir("./data/Gaussian/")
  mkdir("./data/Logistic/")
  mkdir("./data/Quartic/")
  mkdir("./data/Gaussian/HDF5/")
  mkdir("./data/Logistic/HDF5/")
  mkdir("./data/Quartic/HDF5/")
end


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

while isfile("./data/$(thermoname)/HDF5/$(filename)gear.hdf5")
  println("The filename typed already exists in the HDF5 folder. Try another one:")
  filename = string(readline(STDIN))
    filename = filename[1:end-1]
end

parameters = YAML.load(open("parametersgear.yaml"))

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
attrs(file)["rho"] = rho

results =  GearIntegration.run(runtime, rho, dt, T, N, Q, thermoname)


time = results[1]
energy = results[2]
kinetic = results[3]
potential = results[4]
temperature = results[5]
invariant = results[6]
atoms = results[7]   ###For the moment, I don't save it in HDF5 file
zeta = results[8]
nu = results[9]
patom = results[10]
qatom = results[11]

file["t"] = time
attrs(file)["runtime"] = time[end]
file["E"] = energy
file["K"] = kinetic
file["U"] = potential
file["T"] = temperature
file["invariant"] = invariant
file["zeta"] = zeta
file["nu"] = nu
file["patom"] = patom
file["qatom"] = qatom

close(file)

println("The simulation was succesfully done. \nOutput in ./data/$(thermoname)/HDF5/$filename.hdf5")
