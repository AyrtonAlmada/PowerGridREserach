# Compatibility loader for interactive workflows. Use a fresh Julia session when
# switching between legacy and refactored implementations.
include(joinpath(@__DIR__, "src", "PowerGridsFunctions3.jl"))
using .PowerGridsFunctions3
