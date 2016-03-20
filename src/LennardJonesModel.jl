module LennardJonesModel

export initialize, computeforces!, measurekineticenergy, Atom, dim

const dim = 3

type Atom{T}
  r::Array{T,1}
  p::Array{T,1}
  f::Array{T,1}
end

Atom(r) = Atom(r,zeros(dim), zeros(dim))

function makelattice(N::Int, L::Float64)
  latticedistance = L/ceil(cbrt(N))  ##Relationship is N/L^3 = (L/latticedistance)^3/L^3   > rho,. It guarantees that the atoms may be put in the cell.
  println("latticedistance = $latticedistance")
  atoms = Array{Atom,1}(N)
  i = 0
  j = 0
  k = 0
  #First point on the lattice
  latticex = i*latticedistance - 0.5*L
  latticey = j*latticedistance - 0.5*L
  latticez = k*latticedistance - 0.5*L

  atoms[1] = Atom([latticex, latticey, latticez])

  #Lattice point where the atoms are put
  for loop in 2:N
    i+=1
    latticex = i*latticedistance - 0.5*L

    if latticex >= 0.5*L
      i = 0 #Put back to the first position but now second column
      latticex =  i*latticedistance - 0.5*L
      j += 1

    end

    latticey =  j*latticedistance - 0.5*L

    if latticey >= 0.5*L
      j = 0 #Put back to the first position but now third column
      latticey = j*latticedistance - 0.5*L
      k  += 1

    end

    latticez = k*latticedistance - 0.5*L

    atoms[loop] = Atom([latticex, latticey, latticez])

  end

  return atoms
end

@doc """function that creates the array of N atoms in a lattice of side L, with the temperature given by T"""->
function initialize(L::Float64, N::Int64, T::Float64, rho::Float64)
  atoms = makelattice(N, L)
  K = 0 #Kinetic energy

  ##The approach for velocities is adapted from the book "Understanding Molecular Dynamics" (Daan Frenkel)

  sumv = zeros(dim)
  for i=1:N
    for j=1:dim
      v = 2*rand() - 1.0  ##Initial velocity in the interval [-1.0, 1.0]
      sumv[j] += v
      atoms[i].p[j] = v
    end
  end

  ### All the components of the average momentum (p of center of mass) have to be equal to zero. Moreover, the scale factor
  ## guarantees that K = nkT/2
  vaverage =  sumv/(N)
  sumv2new = 0.0

  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] = (atoms[i].p[j] - vaverage[j])
      sumv2new += atoms[i].p[j]^2
    end
  end

  scale = sqrt(dim*(N-1)*T/sumv2new)
  println("scale= $scale")

  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] = scale*atoms[i].p[j]
      K += atoms[i].p[j]^2.
    end
  end

  U = computeforces!(atoms, L)

  #Instantaneous kinetic temperature and energy
  T = K/(dim*(N-1))
  K = K/2

  return atoms, T, K, U
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


@doc """ Determine the interaction force for each pair of particles (i, j)"""->
function computeforces!(atoms::Array{Atom,1}, L::Float64)

  rc = 2.5 #Inner cutoff radius
  N = length(atoms)
  U = 0.0

  for i in 1:N
    for j in 1:dim
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

        r = sqrt(r2)   ####Check how to get rid of this calculation (how to avoid taking the square root)
        r2inverse = 1/r2
        r6inverse = r2inverse * r2inverse * r2inverse


        Vij = 4*r6inverse*(r6inverse - 1.) - 4*(1/rc^12 - 1/rc^6)-(-48./rc^13 + 24./rc^7)*(r-rc)
        fij = 48*r2inverse*r6inverse*(r6inverse - 0.5)  + (-48/rc^13 + 24/rc^7)/r

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
  return U
end


function measurekineticenergy(atoms::Array{Atom,1})
  K = 0.
  N = length(atoms)
  for i in 1:N
    for j in 1:dim
      K +=  atoms[i].p[j]^2.
    end
  end

  K = K/2.
  return K
end

end
