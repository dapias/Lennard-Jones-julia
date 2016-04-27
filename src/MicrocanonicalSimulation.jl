include("./LennardJonesModel.jl")


module MicrocanonicalSimulation

using LennardJonesModel


export verlet!, run

@doc """Function that implements the velocity Verlet algorithm to integrate the equations of motion with step size `dt`"""->
function verlet!(atoms::Array{Atom,1}, dt::Float64, L::Float64)
  N = length(atoms)


  # half-force step
  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] +=  0.5*dt*atoms[i].f[j]
    end
  end

  # full free motion step
  for i in 1:N
    for j in 1:dim
      atoms[i].r[j] += dt*atoms[i].p[j]
    end
  end



  # positions were changed, so recompute the forces
  U = computeenergyandforces!(atoms, L)

  # final force half-step
  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] +=  0.5*dt*atoms[i].f[j]
    end
  end



  return U
end

function initializearrays(numsteps::Int)
  time = Array(Float64, numsteps+1)
  energy = Array(Float64, numsteps+1)
  kinetic = Array(Float64, numsteps+1)
  potential = Array(Float64, numsteps+1)
  temperature = Array(Float64, numsteps+1)

  return time, energy, kinetic, potential, temperature
end




function run(runtime::Float64, rho::Float64, dt::Float64, T::Float64, N::Int64)
  L = cbrt(N/rho)
  numsteps = round(Int, ceil(runtime/dt))
  #Put initial conditions
  atoms, Tinst, K , U = initialize(N, T, rho)
  H = U + K

  time, energy, kinetic, potential, temperature = initializearrays(numsteps)


  println("time")
  println("0.0")

  time[1] = 0.0
  energy[1] = H
  potential[1] = U
  kinetic[1] = K
  temperature[1] = Tinst

  i = 1
  #Perform time steps
  try
    for count in 1:numsteps
      U = verlet!(atoms, dt, L)
      K = measurekineticenergy(atoms)
      H = U + K

      T = K*2/(dim*(N-1))



      time[count+1] = count*dt
      energy[count+1] = H
      potential[count+1] = U
      kinetic[count+1] = K
      temperature[count+1] = T

      #Report results
      #println("$(count*dt), $H, $(U),  $(K), $T")
      println("$(count*dt)")
      i += 1
    end

  catch y   ##If the simulation is stopped, the arrays calculated until this moment are returned
    if isa(y, InterruptException)
      time = time[1:i]
      energy = energy[1:i]
      kinetic = kinetic[1:i]
      potential = potential[1:i]
      temperature = temperature[1:i]
      return time, energy, kinetic, potential, temperature, atoms
    end
  end


  return time, energy, kinetic, potential, temperature, atoms
end



end
