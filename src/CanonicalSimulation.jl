include("./MicrocanonicalSimulation.jl")
include("./ThermostatModel.jl")

module CanonicalSimulation

export run

using LennardJonesModel
using MicrocanonicalSimulation
using ThermostatModel




function thermostatstep!(atoms::Array{Atom,1}, thermo::Thermostat, dt::Float64, N::Int64)

  momentasquare = 0.0

  for i in 1:N
    for j in 1:dim
      momentasquare += (atoms[i].p[j])^2.
    end
  end

  thermo.zeta  += dt/4.*(momentasquare - dim*(N-1)/thermo.beta) ##Taking into account the degrees of freedom

  thermo.nu -= dt/2.*(friction(thermo)/thermo.beta)


  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] = atoms[i].p[j]*exp(dt/2.*(friction(thermo)/thermo.beta))
    end
  end

  thermo.zeta  += dt/4.*(exp(dt*friction(thermo)/thermo.beta)*momentasquare - dim*(N-1)/thermo.beta) ##Taking into account the degrees of freedom

end

function run(runtime::Float64, rho::Float64, dt::Float64, T::Float64, N::Int64, Q::Float64, thermotype::ASCIIString)
  L = cbrt(N/rho)
  numsteps = round(Int, ceil(runtime/dt))
  atoms, Tinst, K , U = initialize(L, N, T, rho)
  H = U + K

  thermomodel = eval(parse(thermotype))

  thermo = thermomodel(Q, 1/T)

  time = Array(Float64, numsteps+1)
  energy = Array(Float64, numsteps+1)
  kinetic = Array(Float64, numsteps+1)
  potential = Array(Float64, numsteps+1)
  temperature = Array(Float64, numsteps+1)
  invariant = Array(Float64, numsteps+1)
  vrandomatom =  Array(Float64, numsteps+1)


  ## Thermo variables
  zetas = Array(Float64, numsteps+1)
  nus =  Array(Float64, numsteps+1)

  ## Selecting an atom to study the distribution of a component of the velocity with time
  randomatom = atoms[2]
  vrandomatom[1] = randomatom.p[2]


  println("time")
  println("0.0")

  time[1] = 0.0
  energy[1] = H
  potential[1] = U
  kinetic[1] = K
  temperature[1] = Tinst
  invariant[1] = H - logrhoextended(thermo)/thermo.beta + thermo.nu*((N-1)*dim)/thermo.beta

  ##Thermo variables
  zetas[1] = thermo.zeta
  nus[1] = thermo.nu

  i = 1
  #Perform time steps
  try
    for count in 1:numsteps

      thermostatstep!(atoms,  thermo, dt, N)
      U =integratestep!(atoms, dt, L)    ##It  can be measured before of applying the thermostat because the thermostat does not affect the positions.
      thermostatstep!(atoms,thermo, dt, N)

      K = measurekineticenergy(atoms)
      H = U + K
      T = K*2/(dim*(N-1))  ##Considering the degrees of freedom



      time[count+1] = count*dt
      energy[count+1] = H
      potential[count+1] = U
      kinetic[count+1] = K
      invariant[count+1] = H - logrhoextended(thermo)/thermo.beta + thermo.nu*((N-1)*dim)/thermo.beta  ##This quantity may have a numerical overflow due
      #to it is equal to extendedrho is equal to exp(-zeta^2). It is better to use the analytical form for log(extended(rho))
      # invariant[count+1] = H + thermo.zeta^2/(2*thermo.Q) + thermo.nu*((N-1)*dim)/thermo.beta
      temperature[count+1] = T

      ####Thermo variables

      zetas[count + 1] = thermo.zeta
      nus[count + 1] = thermo.nu

      vrandomatom[count+1] = randomatom.p[2]


      #############################
      #Report results
      #println("$(count*dt), $H, $(U),  $(K), $T")
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
      vrandomatom = vrandomatom[1:i]

      return time, energy, kinetic, potential, temperature, invariant, atoms, zetas, nus, vrandomatom
    end
  end


  return time, energy, kinetic, potential, temperature, invariant, atoms, zetas,nus, vrandomatom
end

end
