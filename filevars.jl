module PathConstants

    export BASE_DIR, PARTICLE_DIR, SIM_DIR, UTILS_DIR


    const BASE_DIR = raw"/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/filevars.jl"
    const PARTICLE_DIR = joinpath(BASE_DIR, "particle_types")
    const SIM_DIR = joinpath(BASE_DIR, "sims")
    const UTILS_DIR = joinpath(BASE_DIR, "utils")
end