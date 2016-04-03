module ThermostatModel

export Thermostat, Gaussian, Logistic, Quartic, logrhoextended, friction

abstract Thermostat

############Nosé-Hooover ###############3

type Gaussian{T} <: Thermostat
  Q::T ##Parameter that characterizes the distribution
  nu::T
  zeta::T
  beta::T
end

Gaussian(Q, beta) = Gaussian(Q, 0.0, rand(), beta)    ##May change the default 0.0 by rand() for nu
Gaussian(Q, beta, zeta) = Gaussian(Q, 0.0, zeta, beta)


@doc """Logarithm of the distibution in the extended space (f(w))"""->
function logrhoextended(T::Gaussian)
  -T.beta*T.zeta^2/(2.*T.Q)
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Gaussian)
  return -T.beta*T.zeta/T.Q
end
####################################################


############Logistic distribution##############

type Logistic{T} <: Thermostat
  Q::T ##Parameter that characterizes the distribution
  nu::T
  zeta::T
  beta::T
end

Logistic(Q, beta) = Logistic(Q, 0.0, rand(), beta)    ##May change the default 0.0 by rand() for nu
Logistic(Q, beta, zeta) = Logistic(Q, 0.0, zeta, beta)

#With the mean as Q

function logrhoextended(T::Logistic)
  z = T.zeta-T.Q
  distribution = exp(z)/((1+exp(z))^2)
  log(distribution)
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Logistic)
  z = T.zeta-T.Q
  (1 - exp(z))/(1+exp(z))
end

##With the standard deviation as Q (Looks worse than with the mean)

# function logrhoextended(T::Thermostat)
#   z = T.zeta/T.Q
#   distribution = exp(z)/(T.Q*(1+exp(z))^2)
#   log(distribution)
# end

# @doc """Friction coefficient f'(w)/f(w)"""->
# function friction(T::Thermostat)
#   z = T.zeta/T.Q
#   (1 - exp(z))/(T.Q*(1+exp(z)))
# end

##############################################

############Distribution that Fukuda uses###############
type Quartic{T} <: Thermostat
  Q::T ##Parameter that characterizes the distribution
  nu::T
  zeta::T
  beta::T
end

Quartic(Q, beta) = Quartic(Q, 0.0, rand(), beta)    ##May change the default 0.0 by rand() for nu
Quartic(Q, beta, zeta) = Quartic(Q, 0.0, zeta, beta)

# # @doc """Logarithm of the distibution in the extended space (f(w))"""->
function logrhoextended(T::Quartic)
  -T.Q*T.zeta^4
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Quartic)
  return -4*T.Q*T.zeta^3
end
####################################################




###############Cauchy distribution################## (Give a good average but in another T value, the behaviour of peta is also
#monotonous without oscillations. It doesn't seem to work)
##T.Q is the gamma in the distribution (it doesn't seem to work)
# function logrhoextended(T::Thermostat)
#   z = T.zeta^2 + T.Q^2
#   distribution = T.Q/(pi*z)
#   log(distribution)
# end

# @doc """Friction coefficient f'(w)/f(w)"""->
# function friction(T::Thermostat)
#   z = T.zeta^2 + T.Q^2
#   -2*T.zeta/z
# end

##T.Q is x_0 in the distribution
# function logrhoextended(T::Thermostat)
#   z = T.zeta - T.Q
#   distribution = 1/(pi*(1+z^2.))
#   log(distribution)
# end

# @doc """Friction coefficient f'(w)/f(w)"""->
# function friction(T::Thermostat)
#   z = T.zeta - T.Q
#   -2*z/(1+z^2.)
# end
###################################################


end
