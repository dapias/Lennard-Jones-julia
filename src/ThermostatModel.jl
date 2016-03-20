module ThermostatModel

export Thermostat, Gaussian, Logistic, Quartic, logrhoextended, friction

abstract Thermostat

############Nosé-Hooover ###############3

type Gaussian{T} <: Thermostat
  Q::T ##Parameter that characterizes the distribution
  eta::T
  p_eta::T
  beta::T
end

Gaussian(Q, beta) = Gaussian(Q, 0.0, rand(), beta)    ##May change the default 0.0 by rand() for eta
Gaussian(Q, beta, p_eta) = Gaussian(Q, 0.0, p_eta, beta)


@doc """Logarithm of the distibution in the extended space (f(w))"""->
function logrhoextended(T::Gaussian)
  -T.beta*T.p_eta^2/(2.*T.Q)
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Gaussian)
  return -T.beta*T.p_eta/T.Q
end
####################################################


############Logistic distribution##############

type Logistic{T} <: Thermostat
  Q::T ##Parameter that characterizes the distribution
  eta::T
  p_eta::T
  beta::T
end

Logistic(Q, beta) = Logistic(Q, 0.0, rand(), beta)    ##May change the default 0.0 by rand() for eta
Logistic(Q, beta, p_eta) = Logistic(Q, 0.0, p_eta, beta)

#With the mean as Q

function logrhoextended(T::Logistic)
  z = T.p_eta-T.Q
  distribution = exp(z)/((1+exp(z))^2)
  log(distribution)
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Logistic)
  z = T.p_eta-T.Q
  (1 - exp(z))/(1+exp(z))
end

##With the standard deviation as Q (Looks worse than with the mean)

# function logrhoextended(T::Thermostat)
#   z = T.p_eta/T.Q
#   distribution = exp(z)/(T.Q*(1+exp(z))^2)
#   log(distribution)
# end

# @doc """Friction coefficient f'(w)/f(w)"""->
# function friction(T::Thermostat)
#   z = T.p_eta/T.Q
#   (1 - exp(z))/(T.Q*(1+exp(z)))
# end

##############################################

############Distribution that Fukuda uses###############
type Quartic{T} <: Thermostat
  Q::T ##Parameter that characterizes the distribution
  eta::T
  p_eta::T
  beta::T
end

Quartic(Q, beta) = Quartic(Q, 0.0, rand(), beta)    ##May change the default 0.0 by rand() for eta
Quartic(Q, beta, p_eta) = Quartic(Q, 0.0, p_eta, beta)

# # @doc """Logarithm of the distibution in the extended space (f(w))"""->
function logrhoextended(T::Quartic)
  -T.Q*T.p_eta^4
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Quartic)
  return -4*T.Q*T.p_eta^3
end
####################################################




###############Cauchy distribution################## (Give a good average but in another T value, the behaviour of peta is also
#monotonous without oscillations. It doesn't seem to work)
##T.Q is the gamma in the distribution (it doesn't seem to work)
# function logrhoextended(T::Thermostat)
#   z = T.p_eta^2 + T.Q^2
#   distribution = T.Q/(pi*z)
#   log(distribution)
# end

# @doc """Friction coefficient f'(w)/f(w)"""->
# function friction(T::Thermostat)
#   z = T.p_eta^2 + T.Q^2
#   -2*T.p_eta/z
# end

##T.Q is x_0 in the distribution
# function logrhoextended(T::Thermostat)
#   z = T.p_eta - T.Q
#   distribution = 1/(pi*(1+z^2.))
#   log(distribution)
# end

# @doc """Friction coefficient f'(w)/f(w)"""->
# function friction(T::Thermostat)
#   z = T.p_eta - T.Q
#   -2*z/(1+z^2.)
# end
###################################################


end
