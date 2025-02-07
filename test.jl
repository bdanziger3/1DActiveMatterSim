print(1 + 1, "\n")

function f(x, y)
    return x + y
end


print(f(4, 3), "\n")


mutable struct Particle
    mass::Real
    momentofinertia::Real
    x::Real
    v::Real
end

function timestep!(p::Particle, dt::Real, force::Real)
    p.x += (p.v * dt)
    p.v += ((force / p.mass) * dt)
    # return Particle(p.mass, p.momentofinertia, newx, newv)
end

p = Particle(1, 1, 0, 1)

# println(p)

xarray::Array{Real} = []

# xarray.a

println(xarray)

for i in 1:1000
    timestep!(p, .01, -.2)
    push!(xarray, p.x)

    # println(p)
end


include("./plotting/plot_template.jl")

simpleplotxy(0:length(xarray)-1, xarray)


