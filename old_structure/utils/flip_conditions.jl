using Random


"""
Linear approximation of random flipping.

Probability in a timestep dt is dt * T where T is some Temperature parameter related to the reciprocal of some characteristic time.

Valid if dt * T << 1

Returns a function that can be used with a parameter `dt` to simulate random flipping.
"""
function randomlinear(T::Real)::Function
    return (dt::Real -> (rand() < (dt * T)))
end

"""
exponential distribution of polarization flipping events.
Gives a memoryless event distribution

Can be defines with some lambda parameter
or a characteristic time + a probability `p0` of the event happening in a time interval of T:

These two values are related by:

`\\lambda = - \\ln(1 - p0) / T`

pdf function is `pdf(t) = \\lambda * exp(- \\lambda * t)`
cdf function is `pdf(t) = 1 - exp(- \\lambda * t)`
"""
function randomexp(lambda::Real)::Function
    # p0 = (1 - exp(-1 * lambda * dt))
    return (dt::Real -> rand() < (1 - exp(-1 * lambda * dt)))
end

function calclambda(p0::Real, T::Real)::Real
    lambda = - log(1 - p0) / T
    return lambda
end

function randomexp(p0::Real, T::Real)::Function
    lambda = calclambda(p0, T)
    return randomexp(lambda)
end