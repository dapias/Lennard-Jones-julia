include("./src/MicrocanonicalSimulation.jl")

using MicrocanonicalSimulation



#Parameters tested
N = 343
rho = 1.0
T = 0.1
dt = 0.01
runtime = 40.0

time, energy, kinetic, potential, temperature, atoms = MicrocanonicalSimulation.run(runtime, rho, dt, T, N);
