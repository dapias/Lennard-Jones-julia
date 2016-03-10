include("LennardJones.jl")

import LennardJones



#Parameters tested
N = 343
rho = 1.0
T = 0.1
timestep = 0.01
runtime = 40.0

time, energy, kinetic, potential, temperature = LennardJones.run(runtime, rho, timestep, T, N)
