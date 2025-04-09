include("../utils/flip_conditions.jl")

"""
Particle Struct to hold a particle with its spin and velocity.

Particle is currently a point particle with no size or mass.

# Arguments
- `x::Real`: Current location of the particle on a 1D grid. Absolute and non-periodid`x \\in \\Reals`

- `s::Int8`: Initial spin (direction) of particle. S ∈ {1, -1} `S \\in 1`

- `v0::Real`: Initial velocity (in direction S). `v0 \\ge 0`

- `T::Real`: Characteristic time of the particle changing spin.
Probability of changing polarization in time dt is T * dt
"""
mutable struct ParticleState
    x::Real
    S::Int8
    v0::Real
end

mutable struct EnsembleState
    x::Array{Real}
    S::Array{Int8}
    v0::Array{Real}
end

"""
Uses the position parameter x to determine how many periodic cell boundaries the particle
is from its original position.
- `cell::Int`: Counts how many periodic boundaries it has crosses and says which cell it is in
"""
function cell(state::ParticleState, boxwidth::Real)
    celldisplacement = state.x /  boxwidth
    return celldisplacement >= 0 ? ceil(celldisplacement - 0.5) : floor(celldisplacement + 0.5)
end


mutable struct ParticleInteraction
    interactionfunc::Function
end

mutable struct Particle
    name::String
    state::ParticleState
    interaction::ParticleInteraction
end

"""
Definition for a particle that only has its current state (position, speed, polarization)
and a function defining whether it flips polarization.

Flips using a linear condition relying on an input dt and its intrinsic parameter `T`
"""
mutable struct NonInteractingParticle
    state::ParticleState
    flipcondition::Function
    flipconditionparam::Real

    # bulds `flipcondition` to be a function of some `dt` that uses linear conditions based on the parameter T
    NonInteractingParticle(x, S, v0, T) = new(ParticleState(x, S, v0), randomlinear(T), T)
end

"""
function that advances a particle a timestep dt.

Particle moves a distance of `S * v_0 * dt` and then calculates randomly if its polarization changes for the next step 
"""
function timestep!(p::NonInteractingParticle, dt::Real)
    state = p.state
    state.x += (state.S * state.v0 * dt)
    if (p.flipcondition(dt))
        state.S *= -1
    end

end

"""
Definition for a particle that only has its current state (position, speed, polarization)
and a function defining whether it flips polarization.

Flips using a linear condition relying on an input dt and its intrinsic parameter `T`
"""
mutable struct NonInteractingEnsemble
    nparticles::Int
    state::EnsembleState
    flipcondition::Function
    flipconditionparam::Real

    # bulds `flipcondition` to be a function of some `dt` that uses linear conditions based on the parameter T
    NonInteractingEnsemble(x, S, v0, T) = new(length(x), EnsembleState(x, S, v0), randomlinear(T), T)
end

"""
function that advances a particle a timestep dt.

Particle moves a distance of `S * v_0 * dt` and then calculates randomly if its polarization changes for the next step 
"""
function timestep!(ens::NonInteractingEnsemble, dt::Real)
    state = ens.state

    # update positions in place
    state.x .+= (state.S .* state.v0 * dt)

    # update polarities in place
    flips = BitArray(undef, length(state.S))
    arrayflip = x -> ens.flipcondition(dt)
    flips = arrayflip.(state.S)

    state.S[flips .== true] .*= -1

end


"""
Definition for a particle that only has its current state (position, speed, polarization)
and a function defining whether it flips polarization.

Flips using a linear condition relying on an input dt and its intrinsic parameter `T`
"""
mutable struct InteractingParticle
    state::ParticleState
    flipcondition::Function
    flipconditionparam::Real

    # bulds `flipcondition` to be a function of some `dt` that uses linear conditions based on the parameter T
    NonInteractingParticle(x, S, v0, T) = new(ParticleState(x, S, v0), randomlinear(T), T)
end




# function timestep!(p::NonInteractingParticle, dt::Real, boxwidth::Real)
#     state = p.state
#     state.x += (state.S * state.v0 * dt)

#     # check boundary modulus
#     while state.x > boxwidth
#         state.x -= 2 * boxwidth
#     end

#     while state.x < -boxwidth
#         state.x += 2 * boxwidth
#     end


#     if (p.flipcondition(dt))
#         state.S *= -1
#     end

# end
