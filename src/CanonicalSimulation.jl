include("./MicrocanonicalSimulation.jl")
include("./ThermostatModel.jl")

module CanonicalSimulation

export run

using LennardJonesModel
using MicrocanonicalSimulation
using ThermostatModel

"""
Function that corresponds to the operator L_{th} (see equation ___ in paper)
"""
function thermostatstep!(atoms::Array{Atom,1}, thermo::Thermostat, dt::Float64, N::Int64)

  n = dim*(N-1)
  momentasquare = 0.0

  for i in 1:N
    for j in 1:dim
      momentasquare += (atoms[i].p[j])^2.
    end
  end

  thermo.zeta  += dt/4.*(momentasquare - n/thermo.beta)
  thermo.nu -= dt/2.*(friction(thermo)/thermo.beta)

  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] = atoms[i].p[j]*exp(dt/2.*(friction(thermo)/thermo.beta))
    end
  end
  thermo.zeta  += dt/4.*(exp(dt*friction(thermo)/thermo.beta)*momentasquare - n/thermo.beta)
end

"""
Auxiliar function that create the arrays where the data of the simulation will be stored
"""
function initializearrays(numsteps::Int)
  time = Array{Float64}(numsteps+1)
  energy = Array{Float64}(numsteps+1)
  kinetic = Array{Float64}(numsteps+1)
  potential = Array{Float64}(numsteps+1)
  temperature = Array{Float64}(numsteps+1)
  invariant = Array{Float64}(numsteps+1)
  pparticularatom =  Array{Float64}(numsteps+1)
  qparticularatom = Array{Float64}(numsteps+1)

  return time, energy, kinetic, potential, temperature, invariant, pparticularatom, qparticularatom
end

"""
Main function to simulate the system with a given density `rho`, number of particles `N` and equilibrium temperature `T`
during a time `runtime` with a timestep `dt` and a thermostat charactherized by being of a type `thermotype` with parameter
`Q`.
"""
function run(runtime::Float64, rho::Float64, dt::Float64, T::Float64, N::Int64, Q::Float64, thermotype::String)

  n = dim*(N-1)
  numsteps = round(Int, ceil(runtime/dt))
  atoms, Tinst, K , U, L = initialize(N, T, rho)
  H = U + K

  thermomodel = eval(parse(thermotype))
  thermo = thermomodel(Q, 1/T)

  time, energy, kinetic, potential, temperature, invariant, pparticularatom, qparticularatom = initializearrays(numsteps)


  ## Thermo variables
  zetas = Array{Float64}(numsteps+1)
  nus =  Array{Float64}(numsteps+1)

  ## Selecting an atom to study the distribution of a component (y-component) of the velocity and of the position with time
  particularatom = atoms[2]
  pparticularatom[1] = particularatom.p[2]
  qparticularatom[1] = particularatom.r[2]


  println("time")
  println("0.0")

  time[1] = 0.0
  energy[1] = H
  potential[1] = U
  kinetic[1] = K
  temperature[1] = Tinst
  invariant[1] = H - logrhoextended(thermo)/thermo.beta + thermo.nu*n/thermo.beta

  ##Thermo variables
  zetas[1] = thermo.zeta
  nus[1] = thermo.nu

  i = 1
  #Perform time steps
  try
    for count in 1:numsteps

      thermostatstep!(atoms,  thermo, dt, N)
      U =verlet!(atoms, dt, L)    ##It  can be measured before of applying the thermostat because the thermostat does not affect the positions.
      thermostatstep!(atoms,thermo, dt, N)

      K = measurekineticenergy(atoms)
      H = U + K
      T = K*2/n  ##Considering the degrees of freedom



      time[count+1] = count*dt
      energy[count+1] = H
      potential[count+1] = U
      kinetic[count+1] = K
      invariant[count+1] = H - logrhoextended(thermo)/thermo.beta + thermo.nu*n/thermo.beta
      temperature[count+1] = T

      ####Thermo variables

      zetas[count + 1] = thermo.zeta
      nus[count + 1] = thermo.nu

      ####Atom variables
      pparticularatom[count+1] = particularatom.p[2]
      qparticularatom[count+1] = particularatom.r[2]


      ####Report results
      println("$(count*dt)")
      i += 1
    end

  catch y   ##If the simulation is stopped, the arrays calculated until this moment are returende
    if isa(y, InterruptException)
      time = time[1:i]
      energy = energy[1:i]
      kinetic = kinetic[1:i]
      potential = potential[1:i]
      temperature = temperature[1:i]
      invariant = invariant[1:i]
      zetas = zetas[1:i]
      nus = nus[1:i]
      pparticularatom = pparticularatom[1:i]
      qparticularatom =  qparticularatom[1:i]


      return time, energy, kinetic, potential, temperature, invariant, atoms, zetas, nus, pparticularatom, qparticularatom
    end
  end


  return time, energy, kinetic, potential, temperature, invariant, atoms, zetas,nus, pparticularatom, qparticularatom
end

end
