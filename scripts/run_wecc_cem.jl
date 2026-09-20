using TOML, LinearAlgebra, DataFrames, CSV
root=normpath(joinpath(@__DIR__,".."))
include(joinpath(root,"PowerGridsFunctions3.jl"))
cfgpath=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(root,"config/wecc240.toml")
length(ARGS)>=2 || error("Usage: julia --project=. scripts/run_wecc_cem.jl config/wecc240.toml NEW_OUTPUT_DIRECTORY")
out=abspath(ARGS[2]);ispath(out) && error("Output exists; select a new run directory.")
cfg=TOML.parsefile(cfgpath)
BLAS.set_num_threads(cfg["blas_threads"])
branch=joinpath(root,"data/wecc240/processed/branch240E3.csv")
bus=joinpath(root,"data/wecc240/processed/bus240E3.csv")
setup_start=time_ns()
case=load_case(branch,bus)
score=make_surrogate_score(case.df,case.df2;alpha=cfg["alpha"],sigma=cfg["sigma"],
    injection_scale=cfg["injection_scale"],T1=cfg["T1"],T3=cfg["T3"],
    saveat=cfg["saveat"],dt=cfg["internal_dt"],backend=Symbol(cfg["backend"]),
    noise_structure=Symbol(cfg["noise_structure"]))
setup_seconds=(time_ns()-setup_start)/1e9
result=run_cem(score,nrow(case.df);N_final=cfg["final_samples"],
    N_adapt_iter=cfg["samples_per_adaptation_iteration"],K=cfg["adaptation_iterations"],
    rho=cfg["elite_fraction"],kappa=cfg["old_parameter_retention"],
    r0=cfg["proposal_initial_rate"],lambda=cfg["nominal_duration_rate"],gammas=cfg["thresholds"],
    seed_adapt=cfg["seed_adaptation"],seed_estimate=cfg["seed_estimation"])
cfg["input_loading_and_nominal_preparation_s"]=setup_seconds
save_run(out,result;metadata=cfg,input_paths=[branch,bus,cfgpath],
    extra_tables=Dict("bus_map.csv"=>case.bus_map,"branch_map.csv"=>case.branch_map))
println("Saved complete run records to ",out)
