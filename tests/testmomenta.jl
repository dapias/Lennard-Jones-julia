include("../src/CanonicalSimulation.jl")

##Module that checks conservation of the three components of the momentum for the center of mass

using LennardJonesModel
using MicrocanonicalSimulation
using CanonicalSimulation
using FactCheck


#Parameters tested
N = 343
rho = 1.0
T = 0.1
dt = 0.01
runtime = 1.0
Q = 1.0

L = cbrt(N/rho)

atoms, T, K, U = LennardJonesModel.initialize(L, N, T, rho)


facts("Conserved Linear Momenta At the Beginning of the simulation ") do

  momenta = Array(Float64,3)

  for j in 1:dim
    sumv = 0.0
    for i in 1:N
      sumv += atoms[i].p[j]
    end
    momenta[j] = sumv
  end

  @fact abs(momenta[1]) < 1.0e-10 --> true
  @fact abs(momenta[2]) < 1.0e-10 --> true
  @fact abs(momenta[3]) < 1.0e-10 --> true
end

facts("Conserved Linear Momenta After 100 steps (Microcanonical ensemble) ") do

  time, energy, kinetic, potential, temperature, atoms = MicrocanonicalSimulation.run(runtime, rho, dt, T, N)

  momenta = Array(Float64,3)

  for j in 1:dim
    sumv = 0.0
    for i in 1:N
      sumv += atoms[i].p[j]
    end
    momenta[j] = sumv
  end

  @fact abs(momenta[1]) < 1.0e-10 --> true
  @fact abs(momenta[2]) < 1.0e-10 --> true
  @fact abs(momenta[3]) < 1.0e-10 --> true
end

facts("Conserved Linear Momenta After 100 steps (Canonical ensemble) ") do

  time, energy, kinetic, potential, temperature, invariant, atoms = CanonicalSimulation.run(runtime, rho, dt, T, N, Q)

  momenta = Array(Float64,3)

  for j in 1:dim
    sumv = 0.0
    for i in 1:N
      sumv += atoms[i].p[j]
    end
    momenta[j] = sumv
  end

  @fact abs(momenta[1]) < 1.0e-10 --> true
  @fact abs(momenta[2]) < 1.0e-10 --> true
  @fact abs(momenta[3]) < 1.0e-10 --> true
end




