#Julia code to perform a 3D Lennard-Jones simulation with periodic boundary conditions and limits [-L/2, L/2]
# Based on the example givn in the course CHM1464H (Topics in statistical mechanics: The
# Foundations of molecular simulations)

using Distributions

const dim = 3

type Atom{T}
  r::Array{T,1}
  p::Array{T,1}
  f::Array{T,1}
end

Atom(r) = Atom(r,zeros(dim), zeros(dim))

function makelattice(latticedistance::Float64, N::Int, L::Float64)
  system = Array{Atom,1}(N)

  i = 0
  j = 0
  k = 0

  #Lattice point where the atoms are put
  for loop in 0:N-1
    i+= 1
    latticex = i*latticedistance - 0.5*L

    if latticex >= 0.5*L
      i = 0 #Put back to the first position but now second column
      latticex = i*latticedistance - 0.5*L
      j += 1

    end

    latticey = j*latticedistance - 0.5*L

    if latticey >= 0.5*L
      j = 0 #Put back to the first position but now third column
      latticey = j*latticedistance - 0.5*L
      k  += 1

    end

    latticez = k*latticedistance - 0.5*L

    system[loop + 1] = Atom([latticex, latticey, latticez])

  end

  return system

end




function initialize(sytem::Array{Atom,1}, L::Float64, N::Int64, T::Float64)
  latticedistance = L/ceil(cbrt(N))  #We round the cubic root of N to Inf and then divide L by it, taken from the example but it is not clear the physical meaning of this choice.
  scale = sqrt(T)
  K = 0 #Kinetic energy

  normal = Normal()  #Normal distribution with mean zero and standard deviation one.

  #Initial momentum
  for i=1:N
    x = rand(normal, dim)
    system[i].p = scale*x
    K += sum(system[i].p .^2)
  end

  #Intantaneous kinetic temperature and energy

  T = K/(3*N)
  K = K/2

 end


@doc """Function dealing with periodic boundary conditions. Increases or decreases a number d
until is bounded by -L/2 <= d < L/2"""->
function makePeriodic(distance::Float64, L::Float64)
  while distance < -0.5*L
    distance += L
  end

  while distance >= 0.5*L
    distance -= L
  end

  return distance
end


function computeforces(system::Array{Atom,1})

  #Initialize energy to zero
  U = 0
end
