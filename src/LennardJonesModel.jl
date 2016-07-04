module LennardJonesModel

export initialize, computeenergyandforces!, measurekineticenergy, Atom, dim, makePeriodic, rc

const dim = 3
const rc = 2.5 ##Cutoff radius LennardJones potential

"""
        Type with attributes position (r), momentum (p) and force(f), higher derivatives of r: aa, aaa, a4; and the attribute correction which stores the correction value of the force in the Predictor-Corrector scheme
        """
type Atom{T}
    r::Array{T,1}
    p::Array{T,1}
    f::Array{T,1}
    aa::Array{T,1}
    aaa::Array{T,1}
    a4::Array{T,1}
    correction::Array{T,1} 
end

Atom(r) = Atom(r,zeros(dim), zeros(dim), zeros(dim),zeros(dim),zeros(dim), zeros(dim))

"""
        Given a number of particles (N) and a density (rho) it arranges the atoms in a cubic simulation cell taking a fcc unit cell
        """
function makelattice(N::Int, rho::Float64)

    L = cbrt(N/rho)

    if L/2. < rc
        error("The given parameters violate the condition ```Length of simulation cell/2 > rc```. Please decrease the density or increase the number of particles.")
    end

    # Find M large enough to fit N atoms on the simulation cell with unit cell fcc
    # There are 4 atoms that generate the simulation cell by translation

    M = 1
    while (4*M^3. < N)
        M += 1
    end

    latticedistance = L / M   #This lattice distance guarantees that the atoms may be arranged in the cell.

    #4 atomic positions (base) in the fcc cell

    xcell = [0.25,0.75, 0.75, 0.25]
    ycell = [0.25,0.75, 0.25, 0.75]
    zcell = [0.25,0.25, 0.75, 0.75]

    atoms = Array{Atom,1}(N)

    n = 1
    for x in 0:M-1
        for y in 0:M-1
            for z in 0:M-1
                for k in 1:4
                    latticex = (x + xcell[k])*latticedistance - L/2.
                    latticey = (y + ycell[k])*latticedistance - L/2.
                    latticez = (z + zcell[k])*latticedistance -L/2.
                    if n <= N
                        atoms[n] = Atom([latticex, latticey, latticez])
                        n += 1
                    end
                end
            end
        end
    end

    return atoms, L
end

"""
        Creates the simulation cell that contains N atoms and a density rho thermalized at temperature T
        """
function initialize(N::Int64, T::Float64, rho::Float64)

    atoms, L = makelattice(N, rho)
    n = dim*(N -1)  #Degrees of freedom
    K = 0           #Kinetic energy

    sumv = zeros(dim)
    for i=1:N
        for j=1:dim
            v = 2*rand() - 1.0  ##Initial velocity in the interval [-1.0, 1.0]
            sumv[j] += v
            atoms[i].p[j] = v
        end
    end

    # All the components of the momentum of the center of mass have to be equal to zero. The scale factor
    # guarantees that K = nkT/2

    vaverage =  sumv/(N)
    sumv2new = 0.0

    for i in 1:N
        for j in 1:dim
            atoms[i].p[j] = (atoms[i].p[j] - vaverage[j])
            sumv2new += atoms[i].p[j]^2
        end
    end

    scale = sqrt(n*T/sumv2new)

    for i in 1:N
        for j in 1:dim
            atoms[i].p[j] = scale*atoms[i].p[j]
            K += atoms[i].p[j]^2.
        end
    end

    U = computeenergyandforces!(atoms, L)

    #Instantaneous kinetic temperature and energy
    T = K/n
    K = K/2

    return atoms, T, K, U, L
end

"""
        Function dealing with periodic boundary conditions. Increases or decreases by L
        until distance is bounded by -L/2 <= distance < L/2
        """
function makePeriodic(distance::Float64, L::Float64)
    while distance < -0.5*L
        distance += L
    end

    while distance >= 0.5*L
        distance -= L
    end

    return distance
end


"""
        Determines the interaction force for each pair of particles (i, j), save it in each particle's array `f` and computes
        the potential energy of the array of atoms.
        """
function computeenergyandforces!(atoms::Array{Atom,1}, L::Float64)

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

                r = sqrt(r2)   #Check how to get rid of this calculation
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



"""
        Given the array of atoms it computes the kinetic energy
        """
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
