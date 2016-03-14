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
  -T.beta*T.p_eta^2/(2.*T.Q)   ##Q is a parameter that characterizes the distribution (standard deviation?)
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Thermostat)
  return -T.beta*T.p_eta/T.Q
end
####################################################3

end
