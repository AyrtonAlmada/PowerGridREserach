###############################################################
# Power Grid Dynamics Simulation Toolkit
#
# Author: Ayrton Pablo Almada Jimenez
#
# Description
# -------------------------------------------------------------
# This file contains all functions required to simulate and
# analyze power grid dynamics under transmission line faults.
#
# The implementation is based on the swing equation formulation
# and represents the grid as a weighted graph.
#
# The toolkit supports:
#   • Analytical dynamic solutions
#   • Stochastic differential equation simulations
#   • Fault modelling and clearance
#   • Overheating and stability indicators
#   • Monte-Carlo experiments
#
# The grid is described by two main datasets:
#
#   df   → transmission line information
#   df2  → bus/generator parameters
#
###############################################################

using LinearAlgebra
using DataFrames
using CSV
using Statistics
using SparseArrays
using Random
using Distributions
using Roots
using Plots
using ProgressBars

###############################################################
# NETWORK STRUCTURE
###############################################################

"""
    Laplacian(df)

Construct the weighted Laplacian matrix of the transmission network.

# Description

The Laplacian matrix encodes the electrical connectivity between buses.

For a network with incidence matrix `B` and line susceptances `β`,
the Laplacian is defined as

    L = Bᵀ diag(β) B

# Expected DataFrame format

| Column | Description |
|------|-------------|
| Lines | line index |
| From | sending bus |
| To | receiving bus |
| Susceptance | line susceptance |
| Nnodes | number of buses |

# Returns

`L :: Matrix`

Weighted Laplacian matrix.

# Importance

This matrix defines **how power flows propagate through the network**
and is fundamental to all swing-equation simulations.
"""
function Laplacian(df)

    k = nrow(df)
    Nnodes = df.Nnodes[1]

    # incidence matrix
    IN = zeros(Int, k, Nnodes)

    setindex!.(Ref(IN), 1, df.Lines, df.From)
    setindex!.(Ref(IN), -1, df.Lines, df.To)

    B = IN

    # weighted Laplacian
    L = B' * (df.Susceptance .* B)

    return L

end


###############################################################
# SYSTEM MATRICES (SWING EQUATION)
###############################################################

"""
    XiMass(df_lines, df_nodes)

Construct the mass matrix and system dynamics matrix of the
linearized swing equation.

# Swing Equation

    M dω/dt + Dω + Lθ = P
    dθ/dt = ω

Rewritten as first-order system

    A₃ ẋ = A₄ x + b

where

    x = [ω, θ]

# Inputs

`df_lines` : transmission network data

`df_nodes` : generator parameters

Required columns

| Column | Meaning |
|------|-------|
| Inertia | generator inertia |
| Damping | damping coefficient |

# Returns

`A3` : mass matrix  
`A4` : system matrix

# Importance

These matrices represent the **core physics of the grid dynamics**.
All analytical and numerical simulations rely on them.
"""
function XiMass(df_lines, df_nodes)

    n = df_lines.Nnodes[1]

    L = Laplacian(df_lines)

    M = diagm(df_nodes.Inertia)
    D = diagm(df_nodes.Damping)

    Id = Matrix{Float64}(I, n, n)

    ###########################################################
    # MASS MATRIX
    ###########################################################

    A3 = zeros(2n, 2n)

    A3[1:n,1:n] = M
    A3[n+1:2n,n+1:2n] = Id

    ###########################################################
    # SYSTEM MATRIX
    ###########################################################

    A4 = zeros(2n, 2n)

    A4[1:n,1:n] = -D
    A4[1:n,n+1:2n] = -L
    A4[n+1:2n,1:n] = Id

    return A3, A4

end



###############################################################
# ANALYTICAL SOLUTION GENERATORS
###############################################################

"""
    make_f1(...)

Construct a time-dependent function representing the analytical
solution of the system during a fault.

This function precomputes all expensive linear algebra operations
so that evaluating the solution for many time points is fast.

# Returns

`f(t)` — system state at time `t`.
"""
function make_f1(Q, Δ::Diagonal, Y0, λ, A32, b, nozev, n, G)

    T = promote_type(eltype(Q), eltype(Y0))

    Qd = Matrix{T}(Q)
    Δv = Vector{T}(diag(Δ))
    Y0v = Vector{T}(Y0)

    Qinv = inv(Qd)

    A32inv_b = A32 \ b

    qY0 = Qinv * Y0v
    qA32 = Qinv * A32inv_b

    λmod = vcat(λ[1:end-nozev], ones(T,nozev))
    invλmod = one(T) ./ λmod

    expλ = similar(Δv)
    tmp  = similar(qY0)
    out  = similar(qY0)

    return function (t::Real)

        tT = float(t)

        @. expλ = exp(tT * Δv)

        @. tmp = expλ * qY0

        out .= Qd * tmp

        @. tmp = (expλ - one(T)) * invλmod

        if nozev > 0
            tmp[end-nozev+1:end] .= tT
        end

        @. tmp *= qA32

        out .+= Qd * tmp

        return copy(out)

    end
end


"""
    make_f2(...)

Construct the analytical solution generator for the **post-fault**
grid configuration.
"""
function make_f2(Q, Δ::Diagonal, Y1, λ, A31, b, n)

    T = promote_type(eltype(Q), eltype(Y1))

    Qd = Matrix{T}(Q)
    Δv = Vector{T}(diag(Δ))

    Qinv = inv(Qd)

    A31inv_b = A31 \ b

    qY = Qinv * Y1
    qA31 = Qinv * A31inv_b

    λmod = vcat(λ[1:end-1], one(T))
    invλmod = one(T) ./ λmod

    expλ = similar(Δv)
    tmp  = similar(qY)
    out  = similar(qY)

    return function(t::Real)

        tT = T(t)

        @. expλ = exp(tT * Δv)

        @. tmp = expλ * qY
        out .= Qd * tmp

        @. tmp = (expλ - one(T)) * invλmod

        tmp[end] = tT

        @. tmp = tmp * qA31

        out .+= Qd * tmp

        return out

    end

end



###############################################################
# ANALYTICAL SOLUTION GENERATORS
###############################################################

"""
    make_f1(...)

Construct a time-dependent function representing the analytical
solution of the system during a fault.

This function precomputes all expensive linear algebra operations
so that evaluating the solution for many time points is fast.

# Returns

`f(t)` — system state at time `t`.
"""
function make_f1(Q, Δ::Diagonal, Y0, λ, A32, b, nozev, n, G)

    T = promote_type(eltype(Q), eltype(Y0))

    Qd = Matrix{T}(Q)
    Δv = Vector{T}(diag(Δ))
    Y0v = Vector{T}(Y0)

    Qinv = inv(Qd)

    A32inv_b = A32 \ b

    qY0 = Qinv * Y0v
    qA32 = Qinv * A32inv_b

    λmod = vcat(λ[1:end-nozev], ones(T,nozev))
    invλmod = one(T) ./ λmod

    expλ = similar(Δv)
    tmp  = similar(qY0)
    out  = similar(qY0)

    return function (t::Real)

        tT = float(t)

        @. expλ = exp(tT * Δv)

        @. tmp = expλ * qY0

        out .= Qd * tmp

        @. tmp = (expλ - one(T)) * invλmod

        if nozev > 0
            tmp[end-nozev+1:end] .= tT
        end

        @. tmp *= qA32

        out .+= Qd * tmp

        return copy(out)

    end
end


"""
    make_f2(...)

Construct the analytical solution generator for the **post-fault**
grid configuration.
"""
function make_f2(Q, Δ::Diagonal, Y1, λ, A31, b, n)

    T = promote_type(eltype(Q), eltype(Y1))

    Qd = Matrix{T}(Q)
    Δv = Vector{T}(diag(Δ))

    Qinv = inv(Qd)

    A31inv_b = A31 \ b

    qY = Qinv * Y1
    qA31 = Qinv * A31inv_b

    λmod = vcat(λ[1:end-1], one(T))
    invλmod = one(T) ./ λmod

    expλ = similar(Δv)
    tmp  = similar(qY)
    out  = similar(qY)

    return function(t::Real)

        tT = T(t)

        @. expλ = exp(tT * Δv)

        @. tmp = expλ * qY
        out .= Qd * tmp

        @. tmp = (expλ - one(T)) * invλmod

        tmp[end] = tT

        @. tmp = tmp * qA31

        out .+= Qd * tmp

        return out

    end

end



###############################################################
# STOCHASTIC GRID SIMULATION
###############################################################

"""
    SimulationDEFSDE2(...)

Simulate power grid dynamics using a stochastic differential
equation model with fault events.

# Noise

Gaussian noise applied to generator frequency states.

# Returns

DataFrame containing simulated phase angles.
"""
function SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi=0.0,value=0.001)

    n = df.Nnodes[1]

    P1 = df2.PowerInjections
    P2 = Df2.PowerInjections

    A3,A4 = XiMass(df,df2)
    A5,A6 = XiMass(Df,Df2)

    function f(du,u,p,t)

        if (t>=T1)&&(t<=T2)
            du .= A6*u + [P2;zeros(n)]
        else
            du .= A4*u + [P1;zeros(n)]
        end

    end

    function noise(du,u,p,t)

        for i in 1:n
            du[i] = noi
        end

    end

    L1 = Laplacian(df)

    Theta0 = pinv(L1) * P1

    u0 = [zeros(n);Theta0]

    tspan = (0.0,T3)

    prob = SDEProblem(SDEFunction(f,noise; mass_matrix=A5),u0,tspan)

    sol = solve(prob,ImplicitEM(),dtmax=value)

    RES = hcat(sol.u...)

    MAT = transpose(RES)

    DFX = DataFrame(MAT[:,n+1:2n],:auto)

    DFX.Time = range(0,T3,length=size(DFX)[1])

    return DFX

end



###############################################################
# OVERHEATING INDICATORS
###############################################################

"""
    OverheatingIndicator(DFX, df, df2)

Compute the overheating indicator for a simulated trajectory.

# Description
The overheating indicator measures how long transmission lines
operate beyond their thermal limits.

It integrates the duration where the phase difference across
each transmission line exceeds the allowed threshold.

# Returns

Scalar overheating indicator value.
"""
function OverheatingIndicator(DFX,df,df2)

    th = df.ThetaMax

    n = df.Nnodes[1]
    N = nrow(df)

    tdf = AnalysisGamma(DFX,df,df2)

    M = Float32.(Matrix(tdf[2:end,2:end]))

    Δ = Float32(tdf[1,3] - tdf[1,2])

    colsum = vec(sum(M,dims=1)) .* Δ

    keep = colsum .> 0f0

    any(keep) || return 0f0

    s = 0f0

    @inbounds for j in eachindex(keep)

        keep[j] || continue

        for x in @view M[:,j]

            v = x * Δ

            v > 0f0 && (s += v)

        end

    end

    return s

end


"""
    OverheatingIndicatorIndv(DFX,df,df2)

Compute overheating indicators for each individual line.
"""
function OverheatingIndicatorIndv(DFX,df,df2)

    tdf = AnalysisGammaIndv(DFX,df,df2)

    delta = values(tdf[1,2:end])[2] - values(tdf[1,2:end])[1]

    Opt = tdf[2:end,1:end]

    Opt[1:end,2:end] = Opt[1:end,2:end] .* delta

    Opt.S = sum.(eachrow(Opt[:, names(Opt,Real)]))

    rename!(Opt,:column => :"i,j")

    rename!(Opt,:S => :Sij)

    return Opt[:,[1,end]]

end


###############################################################
# ANALYSIS UTILITIES
###############################################################

"""
    Analysis(DFX,df,df2)

Analyze a trajectory and detect line constraint violations.

# Returns

DataFrame containing

• first exit time  
• first return time  
• faulty line  
• violation intervals
"""
function Analysis(DFX,df,df2)

    th = df.ThetaMax

    n = df.Nnodes[1]
    N = nrow(df)

    for i in 1:N

        str = "X$(df.From[i]),$(df.To[i])"

        df3 = DataFrame(Symbol(str)=>DFX[:,df.From[i]]-DFX[:,df.To[i]])

        DFX = hcat(DFX,df3)

    end

    DF = DFX[:,n+2:size(DFX)[2]]

    B = Bool(1)

    mst = zeros(size(DF)[1],1)

    for i in 1:N

        h = Bool.(abs.(DF[:,i]) .<= th[i])

        mst = hcat(mst,h)

        B = B .&& h

    end

    tra = findall(mst[:,Not([1])] .== 0)

    TRAC = hcat(getindex.(tra,1),getindex.(tra,2))

    TRAC = TRAC[sortperm(TRAC[:,1]),:]

    if size(TRAC)[1] > 0

        Time1 = DFX.Time[minimum(TRAC[:,1])]

        alpha = DFX.Time .> Time1

        beta = alpha .&& B

        Time2 = any(beta) ? minimum(DFX.Time[beta]) : NaN

    else

        Time1 = NaN
        Time2 = NaN

    end

    Results = DataFrame(
        FirstExitTime = [Time1],
        FirstReturnTime = [Time2]
    )

    return Results

end


"""
    AnalysisGamma(DFX,df,df2)

Return boolean matrix indicating line constraint violations.
"""
function AnalysisGamma(DFX,df,df2)

    th = df.ThetaMax

    n = df.Nnodes[1]
    N = nrow(df)

    for i in 1:N

        str = "X$(df.From[i]),$(df.To[i])"

        df3 = DataFrame(Symbol(str)=>DFX[:,df.From[i]]-DFX[:,df.To[i]])

        DFX = hcat(DFX,df3)

    end

    DF = DFX[:,n+1:size(DFX)[2]]

    TDF = DataFrame([[names(DF)];collect.(eachrow(DF))],[:column;Symbol.(axes(DF,1))])

    tdf = copy(TDF)

    tdf[2:end,2:end] = (abs.(TDF[2:end,2:end]) .> th)

    return tdf

end


###############################################################
# COMPARISON METRICS
###############################################################

"""
    ComparisonPhase(df,df2,T1,T2,T3,alpha)

Compare analytical solution and approximation for phase angles.
"""
function ComparisonPhase(df,df2,T1,T2,T3,alpha)

    G = []

    for i in 1:nrow(df)

        Df = copy(df)
        Df[i,4] = alpha * Df[i,4]

        ADFX = AnalyticalSolution(df,df2,Df,df2,T1,T2,T3,0.06)

        BDFX = AnalyticalSolutionApprox(df,df2,Df,df2,T1,T2,T3,0.06)

        Analytical = stack(ADFX[:,1:(end-1)])
        Approx = stack(BDFX[:,1:(end-1)])

        push!(G, sqrt(sum((Approx.value - Analytical.value).^2)))

    end

    return maximum(G)

end


"""
    ComparisonOH(df,df2,T1,T2,T3,alpha)

Compare overheating indicators between analytical and approximate models.
"""
function ComparisonOH(df,df2,T1,T2,T3,alpha)

    G = []
    H = []

    for i in 1:nrow(df)

        Df = copy(df)
        Df[i,4] = alpha * Df[i,4]

        ADFX = AnalyticalSolution(df,df2,Df,df2,T1,T2,T3,0.06)

        BDFX = AnalyticalSolutionApprox(df,df2,Df,df2,T1,T2,T3,0.06)

        push!(G,OverheatingIndicator(ADFX,df,df2))
        push!(H,OverheatingIndicator(BDFX,df,df2))

    end

    return maximum(G), maximum(H)

end


###############################################################
# MONTE-CARLO ESTIMATION UTILITIES
###############################################################

"""
    EstimateGENOHI(Ntimes,df,df2,T1,T2,T3)

Estimate overheating indicator distribution using Monte-Carlo
simulation of random line failures.
"""
function EstimateGENOHI(Ntimes,df,df2,T1,T2,T3,alpha=2/3)

    i = rand(1:nrow(df))

    Df = copy(df)
    Df[i,4] = alpha * Df[i,4]

    DFX = AnalyticalSolution(df,df2,Df,df2,T1,T2,T3,0.1)

    oi = OverheatingIndicator(DFX,df,df2)

    An = DataFrame(OverheatingIndicator = [oi])

    for j in 1:Ntimes

        i = rand(1:nrow(df))

        Df = copy(df)
        Df[i,4] = alpha * Df[i,4]

        DFX = AnalyticalSolution(df,df2,Df,df2,T1,T2,T3,0.1)

        oi = OverheatingIndicator(DFX,df,df2)

        push!(An,(oi,))

    end

    return An

end


"""
    EstimateGEN(Ntimes,df,df2,T1,T2,T3)

Run Monte-Carlo experiments collecting failure statistics.
"""
function EstimateGEN(Ntimes,df,df2,T1,T2,T3,alpha=2/3)

    i = rand(1:nrow(df))

    Df = copy(df)
    Df[i,4] = alpha * Df[i,4]

    DFX = AnalyticalSolution(df,df2,Df,df2,T1,T2,T3,0.1)

    oi = OverheatingIndicator(DFX,df,df2)

    An = Analysis(DFX,df,df2)

    An = hcat(An,DataFrame(OverheatingIndicator=[oi]))

    for j in 1:Ntimes

        i = rand(1:nrow(df))

        Df = copy(df)
        Df[i,4] = alpha * Df[i,4]

        DFX = AnalyticalSolution(df,df2,Df,df2,T1,T2,T3,0.1)

        A = Analysis(DFX,df,df2)

        oi = OverheatingIndicator(DFX,df,df2)

        A = hcat(A,DataFrame(OverheatingIndicator=[oi]))

        An = vcat(An,A)

    end

    return An

end