#Julia code to perform a 3D Lennard-Jones simulation with periodic boundary conditions and limits [-L/2, L/2]
# Based on the example given in the course CHM1464H (Topics in statistical mechanics: The
# Foundations of molecular simulations)


module ThermostattedLennardJones

export run, Atom, makelattice

const dim = 3

type Atom{T}
  r::Array{T,1}
  p::Array{T,1}
  f::Array{T,1}
end

type Thermostat{T}
  Q::T ##Parameter that characterizes the distribution
  eta::T
  p_eta::T
  beta::T
end

Thermostat(Q, beta) = Thermostat(Q, 0.0, rand(), beta)
Thermostat(Q, beta, p_eta) = Thermostat(Q, 0.0, p_eta, beta)
Atom(r) = Atom(r,zeros(dim), zeros(dim))


###Extended f(w). In Nosé-Hoover, f(w) is a gaussian function
function extendedrho(w::Float64,  T::Thermostat)
    exp(-T.beta*w^2/(2.*T.Q))
end

##Friction = f'(w)/f(w)
function frictionNH(w::Float64, T::Thermostat)
  ###Nosé-Hoover
  return -T.beta*w/T.Q  ##Q is a parameter that characterizes the distribution
end

################

function makelattice(N::Int, L::Float64, rho::Float64)
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
  atoms = makelattice(N, L, rho)

  K = 0 #Kinetic energy

  ##The approach for velocities is adapted from the book "Understanding Molecular Dynamics" (Daan Frenkel)
  sumv = 0.
  sumv2 = 0.
  #Initial momentum
  for i=1:N
    for j=1:dim
      v = 2*rand() - 1.0
      sumv += v
      sumv2 += v*v
      atoms[i].p[j] = v
    end
  end

  ### The point is that the average momentum (p of center of mass) has to be equal to zero. The scale factor
  ## guarantees that the total kinetic energy is equal to dim/2 times N times kT (k = 1)
  vaverage =  sumv/(dim*N)
  scale = sqrt(dim*T*N/(sumv2 - (sumv)^2./(3*N)))
  println("scale= $scale")



  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] = scale*(atoms[i].p[j] - vaverage)
      K += atoms[i].p[j]^2.
    end
  end

  U = computeforces!(atoms, L)
  #Intantaneous kinetic temperature and energy
  T = K/(dim*(N-1))
  K = K/2

  ###Thermostat variables
  eta = 0.0
  p_eta = rand()

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

        r = sqrt(r2)
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



#The integration will be performed using the Verlet method



# function measureK(atoms::Array{Atom,1})
#   K = 0

#   for i in 1:length(atoms)
#     for j in 1:dim
#       K += atoms[i].p[j]^2.
#     end
#   end

#   K = K/2

#   return K
# end

function thermostatstep!(atoms::Array{Atom,1}, thermo::Thermostat, dt::Float64, N::Int64)

  momentasquare = 0.0

  for i in 1:N
    for j in 1:dim
      momentasquare += (atoms[i].p[j])^2.
    end
  end

  thermo.p_eta  += dt/4.*(-momentasquare + dim*N/thermo.beta)

  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] = atoms[i].p[j]*exp(dt/2.*(-1/thermo.beta*frictionNH(thermo.p_eta,thermo)))
    end
  end

  momentasquare = 0.0

  for i in 1:N
    for j in 1:dim
      momentasquare += (atoms[i].p[j])^2.
    end
  end

  thermo.p_eta  += dt/4.*(-momentasquare + dim*N/thermo.beta)

  thermo.eta += dt/2.*(1/thermo.beta*frictionNH(thermo.p_eta, thermo))



end








function integratestep!(atoms::Array{Atom,1}, dt::Float64, L::Float64)
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
      #  atoms[i].r[j] = makePeriodic(atoms[i].r[j],L)

    end
  end



  # positions were changed, so recompute the forces
  U = computeforces!(atoms, L)

  # final force half-step
  K = 0.
  for i in 1:N
    for j in 1:dim
      atoms[i].p[j] +=  0.5*dt*atoms[i].f[j]
      K +=  atoms[i].p[j]^2.
    end
  end

  K = K/2.

  return K, U
end




function run(runtime::Float64, rho::Float64, dt::Float64, T::Float64, N::Int64, Q::Float64)
  L = cbrt(N/rho)
  numsteps = round(Int, ceil(runtime/dt))
  #Put initial conditions  (##Why so many T's)
  atoms, Tinst, K , U = initialize(L, N, T, rho)
  H = U + K


  thermo = Thermostat(Q, 1/T)

  time = Array(Float64, numsteps+1)
  energy = Array(Float64, numsteps+1)
  kinetic = Array(Float64, numsteps+1)
  potential = Array(Float64, numsteps+1)
  temperature = Array(Float64, numsteps+1)
  invariant = Array(Float64, numsteps+1)



  #Report results
 # println("time, H, U, K, T")
 # println("0.0, $H, $U,  $K, $Tinst")

  println("time")
  println("0.0")

  time[1] = 0.0
  energy[1] = H
  #kineticperparticle[1] = K/N
  #potentialperparticle[1] = U/N
  potential[1] = U
  kinetic[1] = K
  temperature[1] = Tinst
  invariant[1] = H - log(extendedrho(thermo.p_eta, thermo))/thermo.beta


  i = 1
  #Perform time steps
  try
    for count in 1:numsteps
      thermostatstep!(atoms,  thermo, dt, N)
      K, U =integratestep!(atoms, dt, L)
      thermostatstep!(atoms,thermo, dt, N)
      H = U + K
      T = K*2/(dim*(N-1))  ##Considering the degrees of freedom



      time[count+1] = count*dt
      energy[count+1] = H
      #kineticperparticle[count+1] = K/N
      #potentialperparticle[count+1] = U/N
      potential[count+1] = U
      kinetic[count+1] = K
      invariant[count+1] = H - log(extendedrho(thermo.p_eta,thermo))/thermo.beta + thermo.eta*(N*dim)/thermo.beta

      temperature[count+1] = T

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
      return time, energy, kinetic, potential, temperature, invariant, atoms
    end
  end


  return time, energy, kinetic, potential, temperature, invariant, atoms
end

end

