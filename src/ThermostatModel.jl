module ThermostatModel

export Thermostat, logrhoextended, friction

type Thermostat{T}
  Q::T ##Parameter that characterizes the distribution
  eta::T
  p_eta::T
  beta::T
end

Thermostat(Q, beta) = Thermostat(Q, 0.0, rand(), beta)    ##May change the default 0.0 by rand() for eta
Thermostat(Q, beta, p_eta) = Thermostat(Q, 0.0, p_eta, beta)


############Nosé-Hooover ###############33333
@doc """Logarithm of the distibution in the extended space (f(w))"""->
function logrhoextended(T::Thermostat)
  -T.beta*T.p_eta^2/(2.*T.Q)
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Thermostat)
  return -T.beta*T.p_eta/T.Q
end
####################################################


############Logistic distribution##############
# function logrhoextended(T::Thermostat)
#   z = T.p_eta-T.Q
#   distribution = exp(z)/((1+exp(z))^2)
#   log(distribution)
# end

# @doc """Friction coefficient f'(w)/f(w)"""->
# function friction(T::Thermostat)
#   z = T.p_eta-T.Q
#   (1 - exp(z))/(1+exp(z))
# end
##############################################


end
