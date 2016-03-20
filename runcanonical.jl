include("./src/CanonicalSimulation.jl")

using CanonicalSimulation



#Parameters tested
# N = 343
# rho = 1.0
# T = 10.0
# dt = 0.001
# runtime = 40.0
# Q = 2.0

N = 32 #(Value used by Cho and Joannopolous)
#N = 343
rho = 0.8
T = 1.5
dt = 0.0005
runtime = 1.0
Q = 2.0
thermotype = "Gaussian"

time, energy, kinetic, potential, temperature, invariant, atoms, peta, etas, vrandomatom = CanonicalSimulation.run(runtime, rho, dt, T, N, Q,
                                                                                                                   thermotype);
