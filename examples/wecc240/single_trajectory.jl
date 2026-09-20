using CSV, DataFrames, Random
include(joinpath(@__DIR__, "..", "..", "PowerGridsFunctions3.jl"))
root=normpath(joinpath(@__DIR__,"..",".."))
case=load_case(joinpath(root,"data/wecc240/processed/branch240E3.csv"),
               joinpath(root,"data/wecc240/processed/bus240E3.csv"))
# Start with a deterministic smoke test. Nonzero sigma requires an explicit backend.
sys=prepare_system(case.df,case.df2;injection_scale=1.02)
c=prepare_contingency(sys,1;alpha=0.516,sigma=0.0)
DFX=simulate_trajectory(c;T1=0.0,T2=0.7,T3=2.5,saveat=0.01,include_frequency=true)
out=joinpath(root,"results/wecc240/example")
ispath(out) && error("Output already exists: $out")
mkpath(out)
CSV.write(joinpath(out,"surrogate_trajectory.csv"),DFX)
CSV.write(joinpath(out,"line_overload.csv"),line_overload_indicators(DFX,case.df))
CSV.write(joinpath(out,"phase_diagnostics.csv"),phase_difference_diagnostics(DFX,case.df))
CSV.write(joinpath(out,"bus_map.csv"),case.bus_map)
CSV.write(joinpath(out,"branch_map.csv"),case.branch_map)
