include("./LennardJonesModel.jl")
include("./ThermostatModel.jl")

import Base.run

module GearIntegration

export run

using LennardJonesModel
using ThermostatModel


"""
Auxiliar function that create the arrays where the data of the simulation will be stored
"""
function initializearrays(numsteps::Int)
  time = Array(Float64, numsteps+1)
  energy = Array(Float64, numsteps+1)
  kinetic = Array(Float64, numsteps+1)
  potential = Array(Float64, numsteps+1)
  temperature = Array(Float64, numsteps+1)
  invariant = Array(Float64, numsteps+1)
  pparticularatom =  Array(Float64, numsteps+1)
  qparticularatom = Array(Float64, numsteps+1)

  return time, energy, kinetic, potential, temperature, invariant, pparticularatom, qparticularatom
end

function predictorgear!(atom::Atom, Δt::Float64)
    Δt2 = Δt^2.0/2.
    Δt3 = Δt2*Δt/3.
    Δt4 = Δt3*Δt/4.0
    Δt5 = Δt4*Δt/5.0

    atom.r +=  Δt*atom.p + Δt2*atom.f + Δt3*atom.aa + Δt4*atom.aaa + Δt5*atom.a4
    atom.p += Δt*atom.f + Δt2*atom.aa + Δt3*atom.aaa + Δt4*atom.a4
    atom.f += Δt*atom.aa + Δt2*atom.aaa + Δt3*atom.a4
    atom.aa += Δt*atom.aaa + Δt2*atom.a4
    atom.aaa += Δt*atom.a4   
end

function predictorgear!(atoms::Array{Atom,1}, Δt::Float64)
    for atom in atoms
        predictorgear!(atom,Δt)
    end
end

function predictorgear!(atoms::Array{Atom,1}, thermo::Thermostat,  Δt::Float64)
    predictorgear!(atoms, Δt)
    #predictorgear!(thermo,Δt)    
    
    return nothing
end

function computeforcegear!(atoms::Array{Atom,1}, thermo::Thermostat, L::Float64)
    
     #computeforcegear!(thermo,atoms)  ##Checar este orden
    
    N = length(atoms)
    U = 0.0

    for i in 1:N
        for j in 1:dim
            atoms[i].correction[j] = -1.*copy(atoms[i].f[j])
            atoms[i].f[j] = 0.
        end
    end

    for i in 1:N-1
        for j in (i+1):N

            deltax = makePeriodic(atoms[i].r[1] - atoms[j].r[1], L)
            deltay = makePeriodic(atoms[i].r[2] - atoms[j].r[2], L)
            deltaz = makePeriodic(atoms[i].r[3] - atoms[j].r[3], L)

            r2 = deltax*deltax + deltay*deltay + deltaz*deltaz

            if r2 < rc*rc

                r = sqrt(r2)   #Check how to get rid of this calculation
                r2inverse = 1/r2
                r6inverse = r2inverse * r2inverse * r2inverse


                Vij = 4.0*r6inverse*(r6inverse - 1.) - 4*(1/rc^12 - 1/rc^6)-(-48./rc^13 + 24./rc^7)*(r-rc)
                fij = 48.0*r2inverse*r6inverse*(r6inverse - 0.5)  + (-48/rc^13 + 24/rc^7)/r

                U += Vij

                atoms[i].f[1] += fij*deltax 
                atoms[i].f[2] += fij*deltay 
                atoms[i].f[3] += fij*deltaz 

                atoms[j].f[1] -= fij*deltax
                atoms[j].f[2] -= fij*deltay 
                atoms[j].f[3] -= fij*deltaz 
            end

            end
        end
        
    for i in 1:N
        atoms[i].f += friction(thermo)/thermo.beta*atoms[i].p
        atoms[i].correction += atoms[i].f
    end

    return U

end

function correctorgear!(atom::Atom,Δt::Float64)
    k0 = 3/16.
    k1 = 251/360.
    k3 = 11/18.
    k4 = 1/6.
    k5 = 1/60.

    atom.r += k0*Δt^2./2*atom.correction
    atom.p +=  k1*Δt/2.*atom.correction
    atom.aa += k3*atom.correction*6./(2.*Δt)
    atom.aaa += k4*atom.correction*24./(2.*Δt^2.)
    atom.a4 += k5*atom.correction*120./(2.*Δt^3.)    
end

function correctorgear!(thermo::Thermostat,Δt::Float64, n::Int64, atoms::Array{Atom,1})
    
    for atom in atoms
        for i in 1:dim
            thermo.zeta += Δt*(atom.p[i]^2.)
        end
    end
    
    thermo.zeta -= Δt*(n/thermo.beta)
    thermo.nu += Δt*(-friction(thermo)/thermo.beta)
    
    
#     k0 = 251/720.
#     k2 = 11/12.
#     k3 = 1/3.
#     k4 = 1/24.

#     thermo.nu += k0*Δt*thermo.correctionnu
#     thermo.anu +=  k2*thermo.correctionnu*2./(Δt)
#     thermo.aanu += k3*thermo.correctionnu*6./(Δt^2.)
#     thermo.aaanu += k4*thermo.correctionnu*24/(Δt^3.)
    
    
    
#     k0 = 19/120.
#     k1 = 3/4.
#     k3 = 1/2.
#     k4 = 1/12.

#     thermo.zeta += k0*Δt^2./2*thermo.correctionzeta
#     thermo.zetadot +=  k1*Δt/2.*thermo.correctionzeta
#     thermo.aazeta += k3*thermo.correctionzeta*6./(2.*Δt)
#     thermo.aaazeta += k4*thermo.correctionzeta*24/(2.*Δt^2.)
    
    
end

function correctorgear!(atoms::Array{Atom,1}, Δt::Float64)
    for atom in atoms
        correctorgear!(atom,Δt)
    end
end

function correctorgear!(atoms::Array{Atom,1}, thermo::Thermostat, Δt::Float64, n::Int64)
    correctorgear!(atoms,Δt)
    correctorgear!(thermo,Δt, n, atoms)

    return nothing
end

function run(runtime::Float64, rho::Float64, dt::Float64, T::Float64, N::Int64, Q::Float64, thermotype::ASCIIString)

  n = dim*(N-1)
  numsteps = round(Int, ceil(runtime/dt))
  atoms, Tinst, K , U, L = initialize(N, T, rho)
  H = U + K

  thermomodel = eval(parse(thermotype))
  thermo = thermomodel(Q, 1/T)

  time, energy, kinetic, potential, temperature, invariant, pparticularatom, qparticularatom = initializearrays(numsteps)


  ## Thermo variables
  zetas = Array(Float64, numsteps+1)
  nus =  Array(Float64, numsteps+1)

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
  #computeforcegear!(atoms,thermo,L)
  #Perform time steps
  for count in 1:numsteps
        
        predictorgear!(atoms, thermo,dt)
        U = computeforcegear!(atoms,thermo, L)
        correctorgear!(atoms,thermo, dt, n)



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


  return time, energy, kinetic, potential, temperature, invariant, atoms, zetas,nus, pparticularatom, qparticularatom

end


end
