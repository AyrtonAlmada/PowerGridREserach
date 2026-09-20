using Test, Random, LinearAlgebra, DataFrames, CSV, Statistics
include(joinpath(@__DIR__,"..","PowerGridsFunctions3.jl"))

@testset "Source inputs and graph" begin
    root=joinpath(@__DIR__,"..")
    c=load_case(joinpath(root,"data/wecc240/processed/branch240E3.csv"),
                 joinpath(root,"data/wecc240/processed/bus240E3.csv"))
    @test nrow(c.df)==248
    @test nrow(c.df2)==243
    @test c.bus_map.bus_id[1]==1001
    @test norm(Laplatian(c.df)*ones(243))<1e-9
    @test length(initial_state(c.df,c.df2))==486
end

# A two-node deterministic example with known equilibrium and one nonzero mode.
df=DataFrame(Lines=[1],From=[1],To=[2],Susceptance=[2.0],F=[0.3],Nnodes=[2])
df2=DataFrame(Damping=[0.5,0.5],Inertia=[1.0,1.0],PowerInjections=[1.0,-1.0])
@testset "Piecewise propagation and intervals" begin
    s=prepare_system(df,df2)
    c=prepare_contingency(s,1;alpha=0.516,sigma=0.0)
    t=simulate_trajectory(c;T2=0.37,T3=1.0,saveat=0.1)
    @test all(diff(t.Time).>0)
    @test count(==(0.37),t.Time)==1
    @test t.Time[end]==1.0
    @test t.x1[1]≈s.x0[3]
    @test t.x2[1]≈s.x0[4]
    beyond=simulate_trajectory(c;T2=2.0,T3=1.0,saveat=0.1)
    @test beyond.Time[end]==1.0
    same=simulate_trajectory(c;T2=0.0,T3=1.0,saveat=0.1)
    @test maximum(abs.(same.x1 .- same.x1[1]))<1e-8
    # Constant violation must accumulate exactly the duration, not N_grid*dt.
    constpath=DataFrame(x1=[1.,1.,1.],x2=[0.,0.,0.],Time=[0.,0.3,1.])
    @test OverheatingIndicator(constpath,df)==1.0
    shifted=copy(constpath);shifted.x1.+=12.;shifted.x2.+=12.
    @test OverheatingIndicator(shifted,df)==1.0
end
@testset "Noise structure and explicit backend" begin
    s=prepare_system(df,df2)
    c=prepare_contingency(s,1;alpha=0.516,sigma=1.5,noise_structure=:phase_difference)
    @test norm(only(c.G)^2)==0
    @test norm(only(c.G)*vcat(zeros(2),ones(2)))==0
    @test_throws ArgumentError simulate_trajectory(c;T3=1.0)
    a=simulate_trajectory(c;T2=0.2,T3=0.4,backend=:stratonovich_heun,
        dt=0.001,saveat=0.01,rng=Xoshiro(1))
    b=simulate_trajectory(c;T2=0.2,T3=0.4,backend=:stratonovich_heun,
        dt=0.001,saveat=0.01,rng=Xoshiro(1))
    @test isequal(a,b)
    @test a.x1[1]==s.x0[3]
    diagonal=prepare_contingency(s,1;alpha=0.516,sigma=1.5,noise_structure=:frequency_endpoints)
    @test sum(tr(g*g) for g in diagonal.G)≈4.5
end
@testset "Weights, rates, CEM and exact sample counts" begin
    logw=log.([1.,2.,3.,4.]);d=diagnostics_from_logweights(logw,[0,1,0,1])
    @test d.phat≈1.5
    @test d.ESS≈100/30
    @test d.max_normalized_weight≈0.4
    @test d.se≈sqrt(sum(abs2,[0.,2.,0.,4.].-1.5)/12)
    mc=diagnostics_from_logweights(zeros(4),[0,1,0,1])
    @test mc.phat==0.5
    @test mc.ESS==4
    @test mc.max_normalized_weight==0.25
    none=diagnostics_from_logweights(zeros(4),falses(4))
    @test ismissing(none.se)
    # Analytic event depends only on duration, so Q=exp(-lambda*threshold).
    score=(rng,tau,e)->tau
    batch=sample_scores(20000,2.0,[1.0],score;rng=Xoshiro(12))
    @test abs(mean(batch.tau)-0.5)<0.02
    res=cem_exponential(2;score,N=250,K=2,r0=1.0,rng=Xoshiro(13))
    @test res.N_adapt==500
    @test all(res.π.>0)
    @test sum(res.π)≈1.0
    q=estimate_probability(50000,0.5,[0.8,0.2],0.5,
           (rng,tau,e)->tau>2.;rng=Xoshiro(14))
    @test abs(q.phat-exp(-1))<6*q.se
    @test q.ESS*q.max_normalized_weight>=1-1e-10
    @test scenario_filename(9,4,102,7)=="ParaEMTDFX9#4lLamb102T207.csv"
end
