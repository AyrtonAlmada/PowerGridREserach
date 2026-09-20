module PowerGridsFunctions3

# Research code: see docs/MODEL_AND_IMPLEMENTATION.md before changing backends.
# Historical source is preserved byte-for-byte in legacy/source/.
# This module does not run experiments, install packages, or change the global RNG
# or BLAS settings when included.

using LinearAlgebra, SparseArrays, Random, Statistics, Dates, SHA, TOML, Printf
using CSV, DataFrames

export load_case, validate_case, Laplatian, XiMass, initial_state,
       prepare_system, prepare_contingency, simulate_trajectory, AnalyticalSolution,
       line_overload_indicators, OverheatingIndicator, OverheatingIndicatorIndv,
       phase_difference_diagnostics, make_surrogate_score, moment_rhs,
       cem_exponential, sample_scores, estimate_probability, importance_estimate_exp,
       diagnostics_from_logweights, estimates_from_batch, save_run, run_cem,
       file_sha256, scenario_filename

const IMPLEMENTATION_VERSION = "wecc-release-candidate-1"

_check(c, msg) = c || throw(ArgumentError(msg))
function _positive(x, name)
    _check(x isa Real && isfinite(x) && x > 0, "$name must be finite and positive.")
    return Float64(x)
end
function _nonnegative(x, name)
    _check(x isa Real && isfinite(x) && x >= 0, "$name must be finite and nonnegative.")
    return Float64(x)
end
file_sha256(path) = open(io -> bytes2hex(sha256(io)), path)

"""
    load_case(branch_path, bus_path)

Load the supplied processed schema. Preserve ALL rows and the original bus IDs.
`br_x` is converted to `Susceptance = 1/br_x`; no line deduplication, transformer
insertion, power rescaling, or parameter estimation is performed here. Source
units remain the user's input units. See the data dictionary before interpreting
`rate_a` as a physical limit. Return `(df, df2, bus_map, branch_map)`.
"""
function load_case(branch_path::AbstractString, bus_path::AbstractString)
    branch = CSV.read(branch_path, DataFrame)
    bus = CSV.read(bus_path, DataFrame)
    for c in (:index, :f_bus, :t_bus, :br_x, :rate_a)
        _check(c in propertynames(branch), "Branch input lacks $c.")
    end
    for c in (:index, :p, :m, :d)
        _check(c in propertynames(bus), "Bus input lacks $c.")
    end
    ids = Int.(bus.index)
    _check(length(unique(ids)) == length(ids), "Bus IDs must be unique.")
    pos = Dict(id => k for (k,id) in enumerate(ids))
    f_ids, t_ids = Int.(branch.f_bus), Int.(branch.t_bus)
    _check(all(haskey(pos, x) for x in vcat(f_ids, t_ids)), "Unknown branch endpoint.")
    x = Float64.(branch.br_x)
    _check(all(isfinite, x) && all(>(0), x), "br_x must be finite and positive.")
    n, E = nrow(bus), nrow(branch)
    f, t = [pos[x] for x in f_ids], [pos[x] for x in t_ids]
    df = DataFrame(Lines=1:E, From=f, To=t, Susceptance=1 ./ x,
                   F=Float64.(branch.rate_a), Nnodes=fill(n,E))
    df.ThetaMax = df.F ./ df.Susceptance
    df2 = DataFrame(Damping=Float64.(bus.d), Inertia=Float64.(bus.m),
                    PowerInjections=Float64.(bus.p))
    if :va in propertynames(bus)
        df2.Theta0 = Float64.(bus.va) # retained metadata; equilibrium is recomputed
    end
    bus_map = DataFrame(state_index=1:n, bus_id=ids)
    branch_map = DataFrame(line_row=1:E, original_index=Int.(branch.index),
                          from_bus=f_ids, to_bus=t_ids, From=f, To=t,
                          state_label=["X$(f[k])#$(t[k])" for k in 1:E])
    validate_case(df, df2)
    return (; df, df2, bus_map, branch_map)
end

"""Check the all-dynamic-bus ODE schema; zero-inertia algebraic buses are not supported."""
function validate_case(df::AbstractDataFrame, df2::AbstractDataFrame)
    for c in (:Lines,:From,:To,:Susceptance,:F,:Nnodes)
        _check(c in propertynames(df), "df lacks $c.")
    end
    for c in (:Damping,:Inertia,:PowerInjections)
        _check(c in propertynames(df2), "df2 lacks $c.")
    end
    n, E = nrow(df2), nrow(df)
    _check(n > 0 && E > 0, "The network cannot be empty.")
    _check(all(==(n), df.Nnodes), "Nnodes must equal nrow(df2).")
    _check(length(unique(df.Lines)) == E, "Line identifiers must be unique.")
    _check(all(x -> x isa Integer && 1 <= x <= n, df.From), "Invalid From indices.")
    _check(all(x -> x isa Integer && 1 <= x <= n, df.To), "Invalid To indices.")
    _check(all(df.From .!= df.To), "Self-loop branches are not supported.")
    for c in (:Susceptance,:F)
        v = df[!,c]
        _check(all(x -> x isa Real && isfinite(x) && x >= 0, v), "$c must be finite and nonnegative.")
    end
    _check(all(x -> isfinite(x) && x > 0, df2.Inertia), "Inertia must be positive; no silent pseudoinverse of a singular mass matrix.")
    _check(all(x -> isfinite(x) && x >= 0, df2.Damping), "Damping must be nonnegative.")
    _check(all(isfinite, df2.PowerInjections), "Nonfinite injections.")
    return nothing
end

"""Weighted Laplacian. The historical public spelling `Laplatian` is retained."""
function Laplatian(df::AbstractDataFrame)
    n = Int(df.Nnodes[1]); L = zeros(Float64,n,n)
    for r in 1:nrow(df)
        i,j,b = Int(df.From[r]),Int(df.To[r]),Float64(df.Susceptance[r])
        L[i,i] += b; L[j,j] += b; L[i,j] -= b; L[j,i] -= b
    end
    return L
end

"""Return the source's mass and Xi matrices for state order [omega; theta]."""
function XiMass(df, df2)
    validate_case(df,df2)
    n = nrow(df2)
    mass = Matrix(Diagonal(vcat(Float64.(df2.Inertia), ones(n))))
    Xi = zeros(2n,2n)
    Xi[1:n,1:n] = -Matrix(Diagonal(Float64.(df2.Damping)))
    Xi[1:n,n+1:2n] = -Laplatian(df)
    Xi[n+1:2n,1:n] = Matrix{Float64}(I,n,n)
    return mass, Xi
end

"""
Recompute a zero-frequency equilibrium from the actual injections. Each connected
component must balance. No slack correction or power redistribution is implicit.
The Moore--Penrose solution fixes one mean-angle gauge per connected component.
"""
function initial_state(df,df2; rtol=1e-8)
    L = Laplatian(df); p = Float64.(df2.PowerInjections)
    theta = pinv(L) * p
    residual = norm(L*theta-p)
    _check(residual <= rtol*max(norm(p),1.0),
           "No balanced pre-fault equilibrium: residual=$residual. Check component injections.")
    return vcat(zeros(length(p)),theta)
end

# Scalar phi function avoids assuming the zero eigenvalues are last in eigen(A).
_phi(lambda,t) = abs(lambda*t) < 1e-8 ?
    t*(1 + (lambda*t)/2 + (lambda*t)^2/6) : expm1(lambda*t)/lambda

struct AffinePropagator
    A::Matrix{Float64}
    b::Vector{Float64}
    Q::Matrix{ComplexF64}
    Qi::Matrix{ComplexF64}
    lambda::Vector{ComplexF64}
    spectral::Bool
end
function AffinePropagator(A,b)
    A = Matrix{Float64}(A); b = Vector{Float64}(b)
    decomp = eigen(complex.(A)); Q = Matrix{ComplexF64}(decomp.vectors)
    condition = cond(Q)
    ok = isfinite(condition) && condition < 1e10
    Qi = ok ? inv(Q) : zeros(ComplexF64,size(Q))
    return AffinePropagator(A,b,Q,Qi,ComplexF64.(decomp.values),ok)
end
function _affine(p::AffinePropagator,x0,t)
    t == 0 && return Float64.(x0)
    if p.spectral
        v = p.Q * (exp.(p.lambda*t).*(p.Qi*x0) +
                   [_phi(l,t) for l in p.lambda].*(p.Qi*p.b))
        _check(norm(imag.(v)) <= 1e-7*max(norm(real.(v)),1.0), "Unexpected complex solution residue.")
        return real.(v)
    end
    # Rare defective/ill-conditioned eigenbasis: robust but more expensive fallback.
    d = length(x0); aug = zeros(d+1,d+1)
    aug[1:d,1:d] = p.A; aug[1:d,end] = p.b
    return (exp(t*aug)*vcat(x0,1.0))[1:d]
end

struct PreparedSystem
    df::DataFrame
    df2::DataFrame
    n::Int
    A::Matrix{Float64}
    b::Vector{Float64}
    x0::Vector{Float64}
    propagator::AffinePropagator
end
function prepare_system(df,df2; injection_scale=1.0)
    _positive(injection_scale,"injection_scale")
    d2=copy(df2); d2[!,:PowerInjections]=Float64.(df2.PowerInjections).*injection_scale
    mass,Xi=XiMass(df,d2); A=mass\Xi
    b=mass\vcat(Float64.(d2.PowerInjections),zeros(nrow(d2)))
    return PreparedSystem(copy(df),d2,nrow(d2),A,b,initial_state(df,d2),AffinePropagator(A,b))
end

struct PreparedContingency
    system::PreparedSystem
    altered::DataFrame
    altered_buses::DataFrame
    A::Matrix{Float64}
    b::Vector{Float64}
    propagator::AffinePropagator
    G::Vector{Matrix{Float64}}
    noise_structure::Symbol
    sigma::Float64
    changed_rows::Vector{Int}
end
function _prepare_altered(system,Df,Df2,sigma,noise_structure)
    validate_case(Df,Df2)
    _check(nrow(Df)==nrow(system.df) && Df.From==system.df.From && Df.To==system.df.To,
           "Altered and nominal branches must have identical order/endpoints.")
    mass,Xi=XiMass(Df,Df2); A=mass\Xi
    b=mass\vcat(Float64.(Df2.PowerInjections),zeros(system.n))
    changed=findall(Df.Susceptance .!= system.df.Susceptance)
    G=Matrix{Float64}[]; n=system.n
    _nonnegative(sigma,"sigma")
    _check(noise_structure in (:phase_difference,:frequency_endpoints), "Unknown noise_structure.")
    if sigma > 0
        _check(length(changed)==1, "Stochastic correction requires exactly one changed branch.")
        r=only(changed); i,j=Int(Df.From[r]),Int(Df.To[r])
        if noise_structure==:phase_difference
            # This is exactly the BLOCK LOCATION of G in the uploaded Ysol1.
            v=zeros(n);v[i]=1;v[j]=-1
            g=zeros(2n,2n)
            g[1:n,n+1:2n]=sigma.*((1 ./ Float64.(Df2.Inertia)).*(v*v'))
            push!(G,g) # g^2=0; no Stratonovich-to-Ito drift correction
        else
            # Explicit alternative discussed in the manuscript, NOT the uploaded G.
            for endpoint in (i,j)
                g=zeros(2n,2n);g[endpoint,endpoint]=sigma;push!(G,g)
            end
        end
    end
    return PreparedContingency(system,copy(Df),copy(Df2),A,b,AffinePropagator(A,b),
                               G,noise_structure,Float64(sigma),changed)
end
function prepare_contingency(system::PreparedSystem,e::Integer;
                             alpha=0.516,sigma=0.0,noise_structure=:phase_difference)
    _check(1 <= e <= nrow(system.df),"Line index out of range.")
    _check(isfinite(alpha) && 0 <= alpha <= 1,"alpha must lie in [0,1].")
    Df=copy(system.df);Df[e,:Susceptance]*=alpha
    return _prepare_altered(system,Df,system.df2,sigma,noise_structure)
end

function _time_grid(T1,T2,T3,saveat)
    _nonnegative(T1,"T1");_check(isfinite(T2) && T2>=T1,"T2 must be finite and >= T1.")
    _check(isfinite(T3) && T3>T1,"T3 must be finite and > T1.")
    _positive(saveat,"saveat")
    t=Float64.(collect(0:saveat:T3)); push!(t,Float64(T3),Float64(T1))
    T2 <= T3 && push!(t,Float64(T2))
    # Remove near-duplicate switch times without merging a genuinely positive interval.
    sort!(t); out=Float64[]
    for x in t
        if isempty(out) || x-out[end] > 32eps(max(abs(x),1.0))
            push!(out,x)
        end
    end
    for x in (Float64(T1),Float64(min(T2,T3)),Float64(T3))
        k=argmin(abs.(out .- x)); abs(out[k]-x)<=32eps(max(abs(x),1.0)) && (out[k]=x)
    end
    return unique(sort(out))
end

"""
    simulate_trajectory(contingency; T1=0, T2=0.5, T3=10, saveat=0.01,
                        dt=0.001, backend, rng, include_frequency=false)

`T2` is absolute reclosure time; open-phase duration is T2-T1. A reclosure
beyond T3 leaves the system in its disturbed configuration over the horizon.

Backends (explicit for nonzero sigma):
* :deterministic -- matrix-exponential propagation, sigma must be zero.
* :source_random_map -- preserves the uploaded stochastic evaluation formula;
  independent redraws at output times mean this is NOT a Stratonovich SDE solver.
* :stratonovich_heun -- NEW pathwise predictor-corrector integrator using the
  same Wiener increment in predictor and corrector. Requires convergence tests
  and recalibration; does not reproduce historical random-map results.

Times, pre-switch state and reclosure state are stored exactly once. Heun `dt`
and output `saveat` are distinct. No stochastic kicks occur at t=0 for Heun.
"""
function simulate_trajectory(c::PreparedContingency;T1=0.0,T2=0.5,T3=10.0,
        saveat=0.01,dt=0.001,backend::Symbol=:unspecified,
        rng::AbstractRNG=Random.default_rng(),include_frequency=false)
    _check(backend in (:unspecified,:deterministic,:source_random_map,:stratonovich_heun),"Unknown backend.")
    if backend==:unspecified
        c.sigma==0 ? (backend=:deterministic) : throw(ArgumentError("Nonzero sigma requires an explicit backend; read MODEL_AND_IMPLEMENTATION.md."))
    end
    _check(!(backend==:deterministic && c.sigma>0),"Nonzero noise cannot use the deterministic backend.")
    _check(!(backend==:source_random_map && c.noise_structure!=:phase_difference),
           "The source random-map formula supports only the uploaded phase-difference G.")
    _positive(dt,"dt")
    times=_time_grid(T1,T2,T3,saveat); n=c.system.n; Y=zeros(length(times),2n)
    xopen=_affine(c.system.propagator,c.system.x0,T1)
    x=copy(c.system.x0); xclose=copy(xopen)
    # Prebuild stochastic-map propagators; never diagonalize again at each time.
    propG=nothing; propH=nothing
    if backend==:source_random_map && c.sigma>0
        G=only(c.G)
        propG=AffinePropagator(c.A,G*c.b)
        propH=AffinePropagator(c.A,zeros(2n))
    end
    previous=0.0
    for (k,t) in enumerate(times)
        if t < T1
            x=_affine(c.system.propagator,c.system.x0,t)
        elseif t <= T2
            if backend==:stratonovich_heun && c.sigma>0
                # Split exactly at opening. Output times include T1 and T2.
                if t==T1
                    x=copy(xopen)
                else
                    a=max(previous,Float64(T1))
                    while a < t
                        h=min(Float64(dt),t-a)
                        dw=sqrt(h).*randn(rng,length(c.G))
                        f=c.A*x+c.b
                        kick=zeros(length(x))
                        for r in eachindex(c.G); kick .+= (c.G[r]*x).*dw[r]; end
                        xp=x+h*f+kick
                        xc=x+(h/2).*(f+c.A*xp+c.b)
                        for r in eachindex(c.G)
                            xc .+= 0.5.*(c.G[r]*(x+xp)).*dw[r]
                        end
                        x=xc; a+=h
                    end
                end
            elseif backend==:source_random_map && c.sigma>0
                u=t-T1; G=only(c.G)
                # Source formula is deliberately visible and labeled, not called exact.
                x=_affine(c.propagator,xopen,u)
                xi1=(u+1)*randn(rng);xi3=((u+1)^3/3)*randn(rng)
                x .+= xi1.*_affine(propH,G*xopen,u)
                x .+= (xi1-xi3).*_affine(propG,zeros(2n),u)
            else
                x=_affine(c.propagator,xopen,t-T1)
            end
            t==T2 && (xclose=copy(x))
        else
            if T2==T1; xclose=copy(xopen); end
            x=_affine(c.system.propagator,xclose,t-T2)
        end
        all(isfinite,x) || throw(ErrorException("Nonfinite trajectory at t=$t; do not count failed solves as safe events."))
        Y[k,:]=x;previous=t
    end
    out=DataFrame(Y[:,n+1:2n],Symbol.("x".*string.(1:n)))
    if include_frequency
        for j in 1:n;out[!,Symbol("omega$j")]=Y[:,j];end
    end
    out.Time=times
    return out
end

"""Compatibility entry point; dimension is local, and original nominal inputs are not mutated."""
function AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,value=0.006,bs=0;
        backend=:unspecified,noise_structure=:phase_difference,dt=0.001,
        rng=Random.default_rng(),include_frequency=false)
    s=prepare_system(df,df2)
    c=_prepare_altered(s,Df,Df2,bs,noise_structure)
    return simulate_trajectory(c;T1,T2,T3,saveat=value,dt,backend,rng,include_frequency)
end

"""
Cumulative overload duration on actual time intervals (left-endpoint rule).
Uses nominal beta and F as in the manuscript, NOT a thermal model. No endpoint
wrapping, angle recentering, or threshold rescaling is applied. Every threshold
comparison is strict. Full faults and unbalanced voltages need a richer flow model.
"""
function line_overload_indicators(DFX,df; monitored=collect(1:nrow(df)))
    times=Float64.(DFX.Time);_check(length(times)>=2,"At least two output times required.")
    h=diff(times);_check(all(>(0),h),"Time must be strictly increasing; no duplicated reclosure sample.")
    _check(all(isfinite,times),"Nonfinite times.")
    out=DataFrame(line_row=Int[],From=Int[],To=Int[],Sij=Float64[])
    for r in monitored
        _check(1<=r<=nrow(df),"Invalid monitored line index.")
        i,j=Int(df.From[r]),Int(df.To[r]);a=Float64.(DFX[!,Symbol("x$i")]);b=Float64.(DFX[!,Symbol("x$j")])
        _check(all(isfinite,a) && all(isfinite,b),"Nonfinite phase values.")
        flag=abs.(Float64(df.Susceptance[r]).*(a[1:end-1]-b[1:end-1])) .> Float64(df.F[r])
        push!(out,(r,i,j,sum(h.*flag)))
    end
    return out
end
OverheatingIndicator(DFX,df,df2=nothing) = sum(line_overload_indicators(DFX,df).Sij)
function OverheatingIndicatorIndv(DFX,df,df2=nothing)
    t=line_overload_indicators(DFX,df)
    return DataFrame(Symbol("i,j")=>["X$(r.From),$(r.To)" for r in eachrow(t)],:Sij=>t.Sij)
end
function phase_difference_diagnostics(DFX,df;threshold=0.5)
    _positive(threshold,"threshold")
    times=Float64.(DFX.Time);h=diff(times)
    _check(length(h)>0 && all(>(0),h),"Time must be strictly increasing.")
    out=DataFrame(line_row=Int[],max_abs_difference=Float64[],fraction_time_above=Float64[])
    for r in 1:nrow(df)
        d=abs.(DFX[!,Symbol("x$(df.From[r])")]-DFX[!,Symbol("x$(df.To[r])")])
        _check(all(isfinite,d),"Nonfinite phase differences.")
        push!(out,(r,maximum(d),sum(h.*(d[1:end-1].>threshold))/sum(h)))
    end
    return out # Time fractions for THIS path, not ensemble exceedance probabilities.
end

"""Deterministic moment RHS for the declared linear SDE; not a fitted calibration objective."""
function moment_rhs(A,b,Gs,mu,Sigma;convention=:stratonovich)
    _check(convention in (:stratonovich,:ito),"Invalid calculus convention.")
    drift=copy(A)
    if convention==:stratonovich
        for G in Gs;drift .+= 0.5.*(G*G);end
    end
    dmu=drift*mu+b;dSigma=drift*Sigma+Sigma*drift'
    for G in Gs;dSigma .+= G*(Sigma+mu*mu')*G';end
    return dmu,dSigma
end

"""Build an RNG-aware score with explicit physical parameters and lazily cached line models."""
function make_surrogate_score(df,df2;alpha=0.516,sigma=0.0,injection_scale=1.0,
        T1=0.0,T3=10.0,saveat=0.01,dt=0.001,backend=:unspecified,
        noise_structure=:phase_difference)
    system=prepare_system(df,df2;injection_scale)
    cache=Dict{Int,PreparedContingency}() # serial API; not thread-safe by design
    return function(rng::AbstractRNG,tau::Real,e::Integer)
        _nonnegative(tau,"duration")
        c=get!(cache,Int(e)) do
            prepare_contingency(system,e;alpha,sigma,noise_structure)
        end
        path=simulate_trajectory(c;T1,T2=T1+tau,T3,saveat,dt,backend,rng)
        return OverheatingIndicator(path,system.df,system.df2)
    end
end

# =============================================================================
# CEM + ordinary importance sampling. Exponential parameters are always RATES.
# Adaptation follows the supplied unweighted-elite rule, not weighted rare-event CE.
# =============================================================================

function _probs(x,E;positive=true)
    p=Float64.(collect(x));_check(length(p)==E,"Probability-vector length mismatch.")
    _check(all(isfinite,p) && all(>=(0),p) && sum(p)>0,"Invalid probabilities.")
    p./=sum(p)
    positive && _check(all(>(0),p),"Proposal needs positive mass on every nominal category.")
    return p
end
function _categorical(rng,p)
    u=rand(rng);c=0.0
    for i in eachindex(p); c+=p[i];u<c && return i;end
    return lastindex(p)
end

"""
    sample_scores(N, rate, probabilities, score; rng)

Call `score(sample_rng, duration, line_index)`. All simulator randomness must use
sample_rng. Seeds are stored as hex strings, avoiding CSV/Excel UInt64 truncation.
Exceptions are propagated, never replaced with a safe score. Final batches are
serial and explicit; threading requires a separate reproducibility design.
"""
function sample_scores(N::Integer,rate,probabilities,score;rng=Random.default_rng())
    _check(N>=2,"At least two samples required.");r=_positive(rate,"rate")
    p=_probs(probabilities,length(probabilities))
    tau=Vector{Float64}(undef,N);es=Vector{Int}(undef,N);S=Vector{Float64}(undef,N)
    seedhex=Vector{String}(undef,N)
    for k in 1:N
        tau[k]=randexp(rng)/r;es[k]=_categorical(rng,p)
        seed=rand(rng,UInt64);seedhex[k]=string(seed,base=16,pad=16)
        s=score(Xoshiro(seed),tau[k],es[k])
        _check(s isa Real && isfinite(s) && s>=0,"Score must be finite and nonnegative (sample $k).")
        S[k]=s
    end
    return DataFrame(sample=1:N,tau=tau,e=es,score=S,driver_seed_hex=seedhex)
end

"""
Unweighted elite-based adaptation. kappa retains OLD parameters; default .9.
Exactly K iterations, no claimed convergence test. `r0` is a rate. Both learned
parameters must be retained by final estimation. No ambiguous positional booleans.
"""
function cem_exponential(E::Integer;score,N=250,K=20,rho=0.1,kappa=0.9,
        r0=1.0,pi0=fill(1.0/E,E),rng=Random.default_rng())
    _check(E>=1 && N>=2 && K>=1,"E, N, K must be positive (N>=2).")
    _check(0<rho<1 && 0<kappa<1,"rho and kappa must be strictly between 0 and 1.")
    r=_positive(r0,"r0");p=_probs(pi0,E);nelite=max(1,round(Int,rho*N))
    history=DataFrame(iteration=Int[],elite_threshold=Float64[],best_score=Float64[],
                      rate=Float64[],max_phi=Float64[],min_phi=Float64[])
    phistory=DataFrame(iteration=Int[],e=Int[],phi=Float64[])
    pilots=DataFrame[]
    for iteration in 1:K
        batch=sample_scores(N,r,p,score;rng)
        batch.iteration=fill(iteration,N);batch.sampling_rate=fill(r,N)
        batch.sampling_phi=[p[e] for e in batch.e]
        # Deterministic tie rule: descending score, then ascending sample index.
        elite=sortperm(1:N;by=k->(-batch.score[k],k))[1:nelite]
        counts=zeros(E);for k in elite;counts[batch.e[k]]+=1;end
        bar=mean(batch.tau[elite]);_positive(bar,"elite mean duration")
        r=kappa*r+(1-kappa)/bar
        p=kappa.*p+(1-kappa).*(counts./nelite);p./=sum(p)
        _check(all(>(0),p),"Categorical underflow violated support.")
        batch.is_elite=falses(N);batch.is_elite[elite].=true;push!(pilots,batch)
        push!(history,(iteration,minimum(batch.score[elite]),maximum(batch.score),r,maximum(p),minimum(p)))
        for e in 1:E;push!(phistory,(iteration,e,p[e]));end
    end
    return (r=r,π=p,history=history,phi_history=phistory,pilot_samples=vcat(pilots...),
            N_adapt=N*K,N_per_iteration=N,K=K,rho=Float64(rho),kappa=Float64(kappa),r0=Float64(r0))
end

function _logweights(tau,es,r,p,lambda,pnom)
    _positive(r,"proposal rate");_positive(lambda,"nominal rate")
    E=length(p);q=_probs(p,E);pn=_probs(pnom,E)
    _check(length(tau)==length(es),"Sample length mismatch.")
    out=Vector{Float64}(undef,length(tau))
    for k in eachindex(tau)
        e=es[k];_check(1<=e<=E,"Invalid line category.");_nonnegative(tau[k],"duration")
        out[k]=log(lambda)-log(r)+(r-lambda)*tau[k]+log(pn[e])-log(q[e])
    end
    _check(all(isfinite,out),"Nonfinite log likelihood ratio.")
    return out
end

"""Ordinary IS diagnostics; no self-normalization, no guessed uncertainty at zero hits."""
function diagnostics_from_logweights(logw,indicators)
    N=length(logw);_check(N>=2 && length(indicators)==N,"Need aligned batches of size >=2.")
    _check(all(isfinite,logw),"Logweights must be finite.")
    _check(all(x->x==0 || x==1,indicators),"Indicators must be binary.")
    hit=indicators .== 1;h=count(identity,hit)
    w=exp.(logw.-maximum(logw)); sw=sum(w)
    ess=sw^2/sum(abs2,w); maxnorm=maximum(w)/sw
    if h==0
        phat=0.0;se=missing;rse=missing;event_ess=0.0
    else
        m=maximum(logw[hit]);y=zeros(N);y[hit]=exp.(logw[hit].-m)
        my=mean(y);sy=sqrt(sum(abs2,y.-my)/(N*(N-1)))
        phat=exp(m+log(my));se=sy==0 ? 0.0 : exp(m+log(sy));rse=sy/my
        event_ess=sum(y)^2/sum(abs2,y)
        _check(isfinite(phat) && isfinite(se),"Unrepresentable weighted estimate; retain logweights and inspect tail.")
        # phat can exceed 1 in ordinary finite-sample IS; never clip it.
    end
    mean_weight=exp(maximum(logw)+log(sw)-log(N))
    _check(ess<=N*(1+1e-10) && ess*maxnorm>=1-1e-10,"Weight diagnostic identity failed.")
    return (N_final=N,phat=phat,se=se,relative_se=rse,ESS=ess,
            max_normalized_weight=maxnorm,event_ESS=event_ess,hits=h,mean_weight=mean_weight)
end
function importance_estimate_exp(taus,es,indicators,r,π,lambda,E=length(π);
                                nominal_probabilities=fill(1.0/E,E))
    _check(E==length(π),"E/proposal mismatch.")
    logw=_logweights(taus,es,r,π,lambda,nominal_probabilities)
    return merge(diagnostics_from_logweights(logw,indicators),(logweights=logw,))
end
function estimate_probability(N,r,π,lambda,indicator;rng=Random.default_rng(),
                              nominal_probabilities=fill(1.0/length(π),length(π)))
    # A new sample draw on every call; caller supplies a fresh RNG stream after fitting.
    score=(srng,tau,e)->indicator(srng,tau,e) ? 1.0 : 0.0
    batch=sample_scores(N,r,π,score;rng)
    logw=_logweights(batch.tau,batch.e,r,π,lambda,nominal_probabilities)
    batch.logweight=logw;batch.event=batch.score.>0
    return merge(diagnostics_from_logweights(logw,batch.event),(samples=batch,))
end

"""Evaluate multiple strict thresholds on ONE saved final batch; no extra solves."""
function estimates_from_batch(batch;r,π,lambda=0.1,gammas=[0.0,0.5,5.0,10.0],
                               nominal_probabilities=fill(1.0/length(π),length(π)))
    logw=_logweights(batch.tau,batch.e,r,π,lambda,nominal_probabilities)
    out=DataFrame()
    for gamma in gammas
        _nonnegative(gamma,"gamma")
        d=diagnostics_from_logweights(logw,batch.score.>gamma)
        row=merge((gamma=Float64(gamma),),d)
        push!(out,row;cols=:union)
    end
    saved=copy(batch);saved.logweight=logw
    return out,saved
end

"""
Run adaptation and a fresh final batch with distinct RNG streams. Timings are
measured, not extrapolated. The source scorer's lazily cached preparations are
included when first needed; no claim of a warmed-up benchmark is made. Saving
files happens outside the timer and is explicitly excluded from computational cost.
"""
function run_cem(score,E;N_final=35000,N_adapt_iter=250,K=20,rho=0.1,kappa=0.9,
        r0=1.0,lambda=0.1,gammas=[0.0,0.5,5.0,10.0],
        seed_adapt=12345678,seed_estimate=12345679)
    _check(seed_adapt!=seed_estimate,"Use distinct adaptation and estimation seeds.")
    t=time_ns();res=cem_exponential(E;score,N=N_adapt_iter,K,rho,kappa,r0,rng=Xoshiro(seed_adapt))
    adapt=(time_ns()-t)/1e9
    t=time_ns();batch=sample_scores(N_final,res.r,res.π,score;rng=Xoshiro(seed_estimate))
    summary,saved=estimates_from_batch(batch;r=res.r,π=res.π,lambda,gammas)
    estimate=(time_ns()-t)/1e9
    summary.method=fill("CEM",nrow(summary));summary.N_adapt=fill(res.N_adapt,nrow(summary))
    summary.N_total=summary.N_final.+res.N_adapt
    summary.adaptation_time_s=fill(adapt,nrow(summary));summary.estimation_time_s=fill(estimate,nrow(summary))
    summary.total_time_s=fill(adapt+estimate,nrow(summary))
    meta=Dict{String,Any}("seed_adapt"=>string(seed_adapt),"seed_estimate"=>string(seed_estimate),
        "rng"=>"Random.Xoshiro", "N_final"=>N_final,"N_adapt_per_iteration"=>N_adapt_iter,
        "K"=>K,"rho"=>rho,"kappa_old_weight"=>kappa,"initial_rate_s_inv"=>r0,
        "nominal_rate_s_inv"=>lambda,"proposal_rate_s_inv"=>res.r,"event_comparison"=>">",
        "nominal_line_law"=>"uniform_fault_conditioned","no_fault_mass"=>0.0,
        "initial_line_law"=>"uniform","finite_all_weight_second_moment_sufficient"=>(res.r<2lambda),
        "shared_batch_across_thresholds"=>true,"pilot_weighting"=>"unweighted_elites",
        "elapsed_adaptation_s"=>adapt,"elapsed_estimation_and_diagnostics_s"=>estimate,
        "serialization_in_timing"=>false,"warmup_performed_by_driver"=>false)
    return (adaptation=res,samples=saved,summary=summary,metadata=meta)
end

"""Save complete auditable outputs in a NEW directory; never silently overwrite a run."""
function save_run(path::AbstractString,result;metadata=Dict{String,Any}(),input_paths=String[],
                  extra_tables=Dict{String,DataFrame}())
    _check(!ispath(path),"Run path exists; select a new directory to preserve results.")
    mkpath(path)
    CSV.write(joinpath(path,"summary.csv"),result.summary)
    CSV.write(joinpath(path,"final_samples.csv"),result.samples)
    CSV.write(joinpath(path,"adaptation_history.csv"),result.adaptation.history)
    CSV.write(joinpath(path,"adaptation_phi.csv"),result.adaptation.phi_history)
    CSV.write(joinpath(path,"pilot_samples.csv"),result.adaptation.pilot_samples)
    CSV.write(joinpath(path,"proposal.csv"),DataFrame(e=eachindex(result.adaptation.π),
               phi=result.adaptation.π,rate=fill(result.adaptation.r,length(result.adaptation.π))))
    base=copy(result.metadata)
    base["implementation_version"]=IMPLEMENTATION_VERSION;base["julia_version"]=string(VERSION)
    base["utc_recorded"]=string(now(UTC));base["os"]=string(Sys.KERNEL)
    base["arch"]=string(Sys.ARCH);base["julia_threads"]=Threads.nthreads()
    base["logical_cpu_threads"]=Sys.CPU_THREADS;base["system_memory_bytes"]=Sys.total_memory()
    base["blas_threads"]=BLAS.get_num_threads();base["blas_config"]=string(BLAS.get_config())
    base["cpu"]=isempty(Sys.cpu_info()) ? "unknown" : Sys.cpu_info()[1].model
    base["source_sha256"]=file_sha256(@__FILE__)
    base["experiment"]=Dict(string(k)=>v for (k,v) in pairs(metadata))
    base["inputs"]=[Dict("path"=>abspath(p),"sha256"=>file_sha256(p)) for p in input_paths]
    project=Base.active_project()
    if project!==nothing && isfile(project)
        cp(project,joinpath(path,"Project.toml"))
        manifest=joinpath(dirname(project),"Manifest.toml")
        isfile(manifest) && cp(manifest,joinpath(path,"Manifest.toml"))
    end
    open(joinpath(path,"run.toml"),"w") do io;TOML.print(io,base);end
    for (name,table) in extra_tables
        _check(basename(name)==name && endswith(name,".csv"),"Extra table name must be a CSV basename.")
        dest=joinpath(path,name);_check(!ispath(dest),"Extra table would overwrite a run file.")
        CSV.write(dest,table)
    end
    paths=filter(isfile,readdir(path;join=true))
    CSV.write(joinpath(path,"checksums.csv"),DataFrame(file=basename.(paths),sha256=file_sha256.(paths)))
    return path
end

"""Legacy filename convention; integer codes avoid floating-point formatting errors."""
function scenario_filename(from::Integer,to::Integer,lambda_code::Integer,duration_code::Integer)
    return "ParaEMTDFX$(from)#$(to)lLamb$(lambda_code)T2$(lpad(string(duration_code),2,'0')).csv"
end

end # module
