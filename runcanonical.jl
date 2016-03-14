include("./src/CanonicalSimulation.jl")

using CanonicalSimulation



#Parameters tested
N = 343
rho = 1.0
T = 0.1
dt = 0.01
runtime = 40.0
Q = 1.0

time, energy, kinetic, potential, temperature, invariant, atoms, peta, etas = CanonicalSimulation.run(runtime, rho, dt, T, N, Q);
