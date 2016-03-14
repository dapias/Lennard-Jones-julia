include("ThermostattedLennardJones.jl")

import ThermostattedLennardJones



#Parameters tested
N = 343
rho = 1.0
T = 0.1
dt = 0.01
runtime = 10.0
Q = 1.0

time, energy, kinetic, potential, temperature, invariant, atoms = ThermostattedLennardJones.run(runtime, rho, dt, T, N, Q);
