#using Pkg
#Pkg.add(["LinearAlgebra", "DataFrames", "CSV", "SciMLBase", "OrdinaryDiffEq", "StochasticDiffEq", "Statistics", "SparseArrays","Plots","Interact"])
#,"DiffEqProblemLibrary"
#Pkg.add("AbstractPlotting")
#using Pkg
#Pkg.add(["WebIO", "Interact", "IJulia","AbstractPlotting"])
#Pkg.update()
using LinearAlgebra;
using DataFrames;
using CSV;
using SciMLBase;
#using OrdinaryDiffEq;
#using DiffEqProblemLibrary;
#using StochasticDiffEq;
using Statistics;
using SparseArrays;
using PoissonRandom;
using Distributions;
using CrossEntropyMethod;
using Random;
using Plots;
using SimpleDiffEq;
using ExponentialAction;
using FastExpm;
using ExponentialUtilities;
using Printf;
#using WebIO;
#using AbstractPlotting;
#using Plots;
#using Interact;
#using GLMakie;
using Roots;
using ProgressBars;
using BenchmarkTools;
#BLAS.set_num_threads(1)
# Laplatian matrix
# This in how we expect df to be like
# df = DataFrame(Lines = (1:3), From = [1,2,3], To = [2,3,1], Susceptance=[1,1,1], F =[0.441,0.1,0.318], Nnodes=fill(3, 3))
# df.ThetaMax=df.F./df.Susceptance
function Laplatian(df)
    k = nrow(df)
    Nnodes = df.Nnodes[1]
    IN=zeros(Int,k,Nnodes)
    setindex!.(Ref(IN), 1, df.Lines, df.From)
    setindex!.(Ref(IN), -1, df.Lines, df.To)
    B=IN
    L=B'*(df.Susceptance.*B)
    #L=10*L
    return L
end

#Xi matrix
# This in how we expect df to be like
# df = DataFrame(Lines = (1:3), From = [1,2,3], To = [2,3,1], Susceptance=[1,1,1], F =[0.441,0.1,0.318], Nnodes=fill(3, 3))
# df.ThetaMax=df.F./df.Susceptance
# This in how we expect df2 to be like
# df2 = DataFrame(Damping = fill(2, 3), Inertia = 1:3, PowerInjections = [1,-2,1])

function XiMass(df,df2)
    n=df.Nnodes[1] #Number of of buses
    L=Laplatian(df) #Laplacian Matrix
    M=diagm(df2.Inertia) #Inertia Matrix
    D=diagm(df2.Damping) #Damping Matrix
    Zero=0*Matrix{Float64}(Matrix(I, n, n))
    Id=Matrix{Float64}(Matrix(I, n, n))
    A3=0*Matrix{Float64}(Matrix(I, 2*n, 2*n))
    #Mass
    A3[1:n,1:n]=M;
    A3[1:n,n+1:2*n]=Zero;
    A3[n+1:2*n,1:n]=Zero;
    A3[n+1:2*n,n+1:2*n]=Id;
    A4=0*Matrix{Float64}(Matrix(I, 2*n, 2*n))
    Dneg=-1*D
    Lneg=-1*L
    #Xi
    A4[1:n,1:n]=Dneg;
    A4[1:n,n+1:2*n]=Lneg;
    A4[n+1:2*n,1:n]=Id;
    A4[n+1:2*n,n+1:2*n]=Zero;
    return A3,A4
end







# FAST version for make_f1
function make_f1(
    Q, Δ::Diagonal, Y0, λ, A32, b,
    nozev, n, G
)
    # --------------------------------------------------
    # 1. Promote everything to a single concrete type
    # --------------------------------------------------
    T = promote_type(eltype(Q), eltype(Δ.diag), eltype(Y0),
                     eltype(b), eltype(λ), eltype(G))

    Qd  = Matrix{T}(Q)
    Δv  = Vector{T}(diag(Δ))
    Y0v = Vector{T}(Y0)
    bv  = Vector{T}(b)
    Gd  = Matrix{T}(G)
    λv  = Vector{T}(λ)

    # --------------------------------------------------
    # 2. One-time linear algebra
    # --------------------------------------------------
    Qinv = pinv(Qd)
    #A32inv_b = A32 \ bv
    A32inv_b = pinv(A32) * bv

    qY0   = Qinv * Y0v
    qA32  = Qinv * A32inv_b
    qGY0  = Qinv * (Gd * Y0v)
    qGA32 = Qinv * (Gd * A32inv_b)

    # --------------------------------------------------
    # 3. Eigenvalue handling
    # --------------------------------------------------
    λmod = vcat(λv[1:(end-nozev)], ones(T, nozev))
    invλmod = one(T) ./ λmod

    # --------------------------------------------------
    # 4. Preallocated buffers
    # --------------------------------------------------
    expλ = similar(Δv)
    tmp  = similar(qY0)
    out  = similar(qY0)
    acc  = similar(qY0)

    # --------------------------------------------------
    # 5. Time-dependent closure
    # --------------------------------------------------
    return function (t::Real)
        tT = float(t)

        # exp(tΔ)
        @inbounds @. expλ = exp(tT * Δv)

        # ===== deterministic term 1 =====
        @inbounds @. tmp = expλ * qY0
        mul!(out, Qd, tmp)

        # ===== shared diagonal operator =====
        @inbounds @. tmp = (expλ - one(T)) * invλmod
        if nozev > 0
            @inbounds tmp[end-nozev+1:end] .= tT
        end

        # ===== deterministic term 2 =====
        @inbounds @. tmp *= qA32
        mul!(acc, Qd, tmp)
        out .+= acc

        # ===== stochastic scalars =====
        ξ1 = rand(Normal(0, tT+1))
        #ξ2 = randn()
        ξ3 = rand(Normal(0, ((tT+1)^3)/3))
        
        #CAUTION
        #sqrt_t     = sqrt(tT)
        #inv_sqrt_t = inv(ifelse(sqrt_t==0.0,1,sqrt_t))
        sqrt_t=1
        inv_sqrt_t=1
        

        #=
        sqrt_t     = sqrt(tT)/2
        inv_sqrt_t = inv(ifelse(sqrt_t==0.0,1,sqrt_t))/2
        sqrt_t=max(sqrt_t,inv_sqrt_t)
        inv_sqrt_t=sqrt_t
        =#
        #CAUTION
        
        # ===== stochastic term 3 =====
        @inbounds @. tmp = expλ * qGY0
        mul!(acc, Qd, tmp)
        out .+= (ξ1 * sqrt_t) .* acc

        # ===== stochastic term 4 =====
        @inbounds @. tmp = (expλ - one(T)) * invλmod
        if nozev > 0
            @inbounds tmp[end-nozev+1:end] .= tT
        end
        @inbounds @. tmp *= qGA32
        mul!(acc, Qd, tmp)
        #out .+= (ξ2 * sqrt_t) .* acc
        out .+= (ξ1 * sqrt_t) .* acc

        # ===== stochastic term 5 =====
        out .-= (ξ3 * inv_sqrt_t) .* acc

        return copy(out)   # caller-safe
    end
end




# FAST version for make_f2

#=
function make_f2(Q, Δ::Diagonal, Y1, λ, A31, b, n::Int, nozev)
    T = promote_type(eltype(Q), eltype(Δ.diag), eltype(Y1), eltype(b), eltype(λ))
    Qd = Matrix{T}(Q)
    Δv = Vector{T}(diag(Δ))
    Y1v = Vector{T}(Y1)
    bv  = Vector{T}(b)
    λv  = Vector{T}(λ)

    Qinv_d = inv(Qd)
    A31inv_b = A31 \ bv

    qy = Qinv_d * Y1v
    qa31 = Qinv_d * A31inv_b

    λmod = vcat(λv[1:end-1], one(T))
    invλmod = one(T) ./ λmod

    expλ = similar(Δv)
    tmp  = similar(qy)
    out  = similar(qy)

    return function(t::Real)
        tT = T(t)
        @. expλ = exp(tT * Δv)

        @. tmp = expλ * qy
        mul!(out, Qd, tmp)

        @. tmp = (expλ - one(T)) * invλmod
        tmp[end] = tT
        @. tmp = tmp * qa31
        out .+= Qd * tmp

        return out
    end
end
=#


function make_f2(Q, Δ::Diagonal, Y1, λ, A31, b, n, nozev)
    T = promote_type(eltype(Q), eltype(Δ.diag), eltype(Y1), eltype(b), eltype(λ))
    Qd = Matrix{T}(Q)
    Δv = Vector{T}(diag(Δ))
    Y1v = Vector{T}(Y1)
    bv  = Vector{T}(b)
    λv  = Vector{T}(λ)

    Qinv = inv(Qd)
    A31inv_b = A31 \ bv

    qy = Qinv * Y1v
    qa31 = Qinv * A31inv_b


    λmod = vcat(λv[1:(end-nozev)], ones(T, nozev))
    invλmod = one(T) ./ λmod

    expλ = similar(Δv)
    tmp  = similar(qy)
    out  = similar(qy)

    return function(t::Real)
        tT = T(t)
        @. expλ = exp(tT * Δv)

        @. tmp = expλ * qy
        mul!(out, Qd, tmp)

        @. tmp = (expλ - one(T)) * invλmod


        if nozev > 0
            @. tmp[end-nozev+1:end] .= tT
        end

        
        #tmp[end] = tT
        @. tmp = tmp * qa31
        out .+= Qd * tmp

        return out
    end
end


function Ysol1(df,df2,Df,Df2,bs=0)
    #During Fault
    n=df.Nnodes[1]; #Number of of buses    
    P1=df2.PowerInjections; #Vector of power injections 1
    L1=Laplatian(df); #Laplacian Matrix F
    A32,A42=XiMass(Df,Df2);
    #Afault=pinv(A32)*A42
    #j=findall(Df[:,4]-df[:,4].!=0)[1]
    if size(findall(Df[:,4]-df[:,4].!=0))[1]==0
        j=1
    else
        j=findall(Df[:,4]-df[:,4].!=0)[1]
    end
    Theta0=pinv(L1)*P1; #Steady-state solution 
    #Theta0=L1 \ P1
    Y0=[zeros(n);Theta0];
    #P1=df2.PowerInjections; #Vector of power injections 1
    b=[P1;zeros(n)];

    XIEig=eigen(pinv(A32)*A42);
    λ = XIEig.values;
    Q =  XIEig.vectors;
    Δ = Diagonal(λ);

    nozev=sum(abs.(λ).<1e-10)


    M=diagm(df2.Inertia)

    G=0*Matrix{Float64}(Matrix(I, 2*n, 2*n))

    G[1:n,n+1:2*n]=bs*pinv(M)*((I[1:n,Df[j,:][2]]-I[1:n,Df[j,:][3]])*transpose(I[1:n,Df[j,:][2]]-I[1:n,Df[j,:][3]]))
    
    #KEEP THIS SAFE
    
    #f1 = (t) -> Q*exp(t*Δ)*pinv(Q)*Y0 + Q*((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-nozev],fill(1,nozev))))+Diagonal(vcat(zeros(2*n-nozev),fill(t,nozev))))*pinv(Q)*(pinv(A32)*b)

    f1=make_f1(Q, Δ, Y0, λ, A32, b, nozev, n, G)
    return f1
end




function Ysol2(df,df2,Df,Df2,Y1)
    #Post Fault
    n=df.Nnodes[1]; #Number of of buses    
    P1=df2.PowerInjections; #Vector of power injections 1
    #A31 Mass, A41 Xi
    A31,A41=XiMass(df,df2);
    #Aclear=pinv(A31)*A41
    b=[P1;zeros(n)];

    XIEig=eigen(pinv(A31)*A41);
    λ = XIEig.values;
    Q =  XIEig.vectors;
    Δ = Diagonal(λ);

    nozev=sum(abs.(λ).<1e-10)

    #KEEP THIS SAFE

    #f2 = (t) -> Q*exp(t*Δ)*pinv(Q)*Y1 + Q*((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-1],1)))+Diagonal(vcat(zeros(2*n-1),t)))*pinv(Q)*(pinv(A31)*b)    
    #f2 = make_f2(Q, Δ, Y1, λ, A31, b, n)
    f2 = make_f2(Q, Δ, Y1, λ, A31, b, n, nozev)

    return f2
end

function AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,value=0.006,bs=0)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    f1=Ysol1(df,df2,Df,Df2,bs)
    Tempo1=[T1:value:T2;]
    Tempo2=[T2:value:T3;]
    RES=hcat(real.(f1.(Tempo1))...)
    Y1=real(f1(T2))
    f2=Ysol2(df,df2,Df,Df2,Y1)
    MAT=hcat(RES,hcat(real.(f2.(Tempo2.-T2))...))
    J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    DFX.Time = range(0, T3, length=size(DFX)[1])
    return DFX
end

function AnalyticalSolution2(df,df2,Df,Df2,T1,T2,T3,value=0.006,bs=0)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    Tempo1 = range(T1, T2; step = value)
    Tempo2 = range(T2, T3; step = value)
    
    nT1 = length(Tempo1)
    nT2 = length(Tempo2)

    f1=Ysol1(df,df2,Df,Df2)
    
    # Problem dimension
    m = length(f1(T1))   # or known size = 2n
    
    # ------------------------------------------------------------
    # Preallocate result matrix
    # Rows = time, Columns = state variables
    # ------------------------------------------------------------
    Y = Matrix{Float64}(undef, nT1 + nT2, m)
    
    # ------------------------------------------------------------
    # Fill block 1 (f1)
    # ------------------------------------------------------------
    @inbounds for (i, t) in enumerate(Tempo1)
        Y[i, :] .= real.(f1(t))
    end
    
    # Save final state
    Y1 = view(Y, nT1, :)
    
    # ------------------------------------------------------------
    # Build f2 once
    # ------------------------------------------------------------
    f2 = Ysol2(df, df2, Df, Df2, Y1)
    
    # ------------------------------------------------------------
    # Fill block 2 (f2)
    # ------------------------------------------------------------
    @inbounds for (i, t) in enumerate(Tempo2)
        Y[nT1 + i, :] .= real.(f2(t - T2))
    end
    
    DFX = DataFrame(Y[:, n+1:2n], :auto)
    DFX.Time = range(0, T3; length = size(DFX, 1))
    return DFX
end


function Ysol1xi(df,df2,Df,Df2)
    #During Fault
    n=df.Nnodes[1]; #Number of of buses    
    P1=df2.PowerInjections; #Vector of power injections 1
    L1=Laplatian(df); #Laplacian Matrix F
    #M1=diagm(df2.Inertia); #Inertia Matrix F
    #D1=diagm(df2.Damping); #Damping Matrix F
    #L2=Laplatian(Df); #Laplacian Matrix G
    #M2=diagm(Df2.Inertia); #Inertia Matrix G
    #D2=diagm(Df2.Damping); #Damping Matrix G
    #A32 Mass, A42 Xi
    A32,A42=XiMass(Df,Df2);

    #First Attempt
    #-------------------------------------------------------------------------------------------------------------------------------------------------------#
    #UPDATE THIS PART TO SAMPLE OVER DIFFRENT BASE SOLUTIONS
    #-------------------------------------------------------------------------------------------------------------------------------------------------------#
    #PI=CSV.read("PowerInjections.csv", DataFrame; delim = ',')
    #PI[!,"x33"].=P1
    #i=rand(1:size(PI)[2],1)[1]
    #i=rand(29:size(PI)[2],1)[1]
    #P1=PI[:,i]
    #-------------------------------------------------------------------------------------------------------------------------------------------------------#
    #UPDATE THIS PART TO SAMPLE OVER DIFFRENT BASE SOLUTIONS
    #-------------------------------------------------------------------------------------------------------------------------------------------------------#


    #Second Attempt
    #-------------------------------------------------------------------------------------------------------------------------------------------------------#
    #UPDATE THIS PART TO SAMPLE OVER DIFFRENT BASE SOLUTIONS
    #-------------------------------------------------------------------------------------------------------------------------------------------------------#
    #=
    b=1:n
    gen=[1, 3, 5, 7, 8, 9, 10, 11, 13, 15, 16, 19, 22, 23, 26, 27, 28, 29, 30, 31, 32]
    loa=setdiff(b,gen)
    Gamm=rand(size(gen)[1]) #Random Particiaption Factors
    Gamm./=sum(Gamm)
    deltPL=0.1*rand(size(loa)[1]) #Random Changes in the load demand
    deltP=sum(deltPL)
    DF1=DataFrame(ind=gen,Gamm=Gamm)
    DF2=DataFrame(ind=loa,deltPL=deltPL)
    DF=sort(outerjoin(DF1,DF2, on = :ind),:ind)
    GENS=ismissing.(DF.deltPL).*b;
    PIstar=P1+[ifelse(i ∈ GENS,DF.Gamm[i]*deltP,-DF.deltPL[i]) for i in 1:32]
    P1=PIstar
    =#
    #-------------------------------------------------------------------------------------------------------------------------------------------------------#
    #UPDATE THIS PART TO SAMPLE OVER DIFFRENT BASE SOLUTIONS
    #-------------------------------------------------------------------------------------------------------------------------------------------------------#

    
    
    Theta0=pinv(L1)*P1; #Steady-state solution
    Y0=[zeros(n);Theta0];
    #P1=df2.PowerInjections; #Vector of power injections 1
    b=[P1;zeros(n)];

    XIEig=eigen(pinv(A32)*A42);
    λ = XIEig.values;
    Q = XIEig.vectors;
    Δ = Diagonal(λ);

    nozev=sum(abs.(λ).<1e-10)

    #one zero eigenvalue

    #f1 = (t) -> Q*exp(t*Δ)*pinv(Q)*Y0 + Q*((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-1],1)))+Diagonal(vcat(zeros(2*n-1),t)))*pinv(Q)*(pinv(A32)*b)

    #more than one zero eigenvalue

    f1 = (t) -> exp(t*Δ)*pinv(Q)*Y0 + ((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-nozev],fill(1,nozev))))+Diagonal(vcat(zeros(2*n-nozev),fill(t,nozev))))*pinv(Q)*(pinv(A32)*b)

    
    
    #f1 = (t) -> exp(t*pinv(A32)*A42)*Y0+(exp(t*pinv(A32)*A42)-I)*pinv(pinv(A32)*A42)*pinv(A32)*b;
    #f1 = (t) -> (fastExpm(t*pinv(A32)*A42))*Y0+(fastExpm(t*pinv(A32)*A42)-I)*pinv(pinv(A32)*A42)*pinv(A32)*b;
    return f1,Q
end

function Ysol2xi(df,df2,Df,Df2,Y1)
    #Post Fault
    n=df.Nnodes[1]; #Number of of buses    
    P1=df2.PowerInjections; #Vector of power injections 1
    #A31 Mass, A41 Xi
    A31,A41=XiMass(df,df2);
    b=[P1;zeros(n)];

    XIEig=eigen(pinv(A31)*A41);
    λ = XIEig.values;
    Q =  XIEig.vectors;
    Δ = Diagonal(λ);

    #nozev=sum(abs.(λ).<1e-10)

    #one zero eigenvalue

    f2 = (t) -> exp(t*Δ)*inv(Q)*Y1 + ((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-1],1)))+Diagonal(vcat(zeros(2*n-1),t)))*inv(Q)*(pinv(A31)*b)    

    #more than one zero eigenvalue

    #f2 = (t) -> Q*exp(t*Δ)*pinv(Q)*Y1 + Q*((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-nozev],fill(1,nozev))))+Diagonal(vcat(zeros(2*n-nozev),fill(t,nozev))))*pinv(Q)*(pinv(A31)*b)

    
    #f2 = (t) -> exp(t*pinv(A31)*A41)*Y1+(exp(t*pinv(A31)*A41)-I)*pinv(pinv(A31)*A41)*pinv(A31)*b;
    #f2 = (t) -> (fastExpm(t*pinv(A31)*A41))*Y1+(fastExpm(t*pinv(A31)*A41)-I)*pinv(pinv(A31)*A41)*pinv(A31)*b;    
    return f2,Q
end

function Ysol0(df,df2,Df,Df2)
    n=df.Nnodes[1]; #Number of of buses    
    P1=df2.PowerInjections; #Vector of power injections 1
    L1=Laplatian(df); #Laplacian Matrix F
    
    Theta0=pinv(L1)*P1; #Steady-state solution
    Y0=[zeros(n);Theta0];
    b=[P1;zeros(n)];
    
    #NoFault
    A31,A41=XiMass(df,df2);
    XIEig1=eigen(pinv(A31)*A41);
    Q1 = XIEig1.vectors;
    λ  = XIEig1.values;
    Δ = Diagonal(λ);
    
    Q=Q1
    IQ=inv(Q)

    #Perturbation
    A32,A42=XiMass(Df,Df2);
    V=pinv(A32)*A42-pinv(A31)*A41
    λper=λ+diag(IQ*V*Q)
    Δper = Diagonal(λper)

    A=reshape(repeat(transpose.(λ), outer = 2*n),(2*n,2*n))
    AT=transpose(A)
    DAT=inv.(AT-A+I)-I
    D=(IQ*V*Q).*DAT
    Qpert=Q+Q*D
    
    f1 = (t) -> Qpert*exp(t*Δper)*pinv(Qpert)*Y0 + Qpert*((exp(t*Δper)-I).*pinv(Diagonal(vcat(λper[1:end-1],1)))+Diagonal(vcat(zeros(2*n-1),t)))*pinv(Qpert)*(pinv(A31)*b)
    
    
    return f1
end


function OverHeatingIndicatorXi(df,df2,Df,Df2,T1,T2,T3)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    #f1=Ysol1(df,df2,Df,Df2)
    #Y1=real(f1(T2))
    xi1,Q1=Ysol1xi(df,df2,Df,Df2)
    Y1=Q1*xi1(T2)
    xi2,Q2=Ysol2xi(df,df2,Df,Df2,Y1)
    
    th=df.ThetaMax
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    OHI=[]
    for i in 1:N
        I=df.From[i]
        J=df.To[i]
        fijA = (t) -> abs(real.(Q1[n+I,:]'*xi1(t)-Q1[n+J,:]'*xi1(t)))-th[i]
        rootsA = find_zeros(fijA, T1, T2)
        fijB = (t) -> abs(real.(Q2[n+I,:]'*xi2(t-T2)-Q2[n+J,:]'*xi2(t-T2)))-th[i]
        rootsB = find_zeros(fijB, T2, T3)
        roots=unique(vcat(rootsA,rootsB))
        #roots=vcat(rootsA,rootsB)
        if size(roots)[1]>=2
            push!(OHI,sum(diff(roots)[1:2:end]))
        else
            push!(OHI,0)
        end
    end
    
    return(sum(OHI)) 
end

function Ysol00(df,df2,Df,Df2,Bf,Bf2)
    n=df.Nnodes[1]; #Number of of buses    
    P1=df2.PowerInjections; #Vector of power injections 1
    L1=Laplatian(df); #Laplacian Matrix F
    
    Theta0=pinv(L1)*P1; #Steady-state solution
    Y0=[zeros(n);Theta0];
    b=[P1;zeros(n)];
    
    #NoFault
    A31,A41=XiMass(df,df2);
    XIEig1=eigen(pinv(A31)*A41);
    Q1 = XIEig1.vectors;
    λ  = XIEig1.values;
    Δ = Diagonal(λ);
    
    Q=Q1
    IQ=inv(Q)

    #Single Phase
    A32,A42=XiMass(Df,Df2);
    V=pinv(A32)*A42-pinv(A31)*A41
    λper=λ+diag(IQ*V*Q)
    Δper = Diagonal(λper)

    A=reshape(repeat(transpose.(λ), outer = 2*n),(2*n,2*n))
    AT=transpose(A)
    DAT=inv.(AT-A+I)-I
    D=(IQ*V*Q).*DAT
    Qpert=Q+Q*D

    #Three Phase
    A33,A43=XiMass(Bf,Bf2);
    V1=pinv(A33)*A43-pinv(A32)*A42
    λper2=λper+diag(pinv(Qpert)*V1*Qpert)
    Δper2 = Diagonal(λper2)

    
    A2=reshape(repeat(transpose.(λper), outer = 2*n),(2*n,2*n))
    AT2=transpose(A2)
    DAT2=inv.(AT2-A2+I)-I
    D=(pinv(Qpert)*V1*Qpert).*DAT2
    Qpert2=Qpert+Qpert*D
    
    f1 = (t) -> Qpert2*exp(t*Δper2)*pinv(Qpert2)*Y0 + Qpert2*((exp(t*Δper2)-I).*pinv(Diagonal(vcat(λper2[1:end-1],1)))+Diagonal(vcat(zeros(2*n-1),t)))*pinv(Qpert2)*(pinv(A31)*b)
    
    
    return f1
end

function Ysol000(df,df2,i,m)
    n=df.Nnodes[1]; #Number of of buses    
    P1=df2.PowerInjections; #Vector of power injections 1
    L1=Laplatian(df); #Laplacian Matrix F
    
    Theta0=pinv(L1)*P1; #Steady-state solution
    Y0=[zeros(n);Theta0];
    b=[P1;zeros(n)];

    #NoFault
    A31,A41=XiMass(df,df2);
    XIEig1=eigen(pinv(A31)*A41);
    Q1 = XIEig1.vectors;
    λ  = XIEig1.values;
    Δ = Diagonal(λ);
    
    Q=Q1
    IQ=pinv(Q)
    
    Df=copy(df)
    Df[i,4]=((m-1)/(m))*Df[i,4]
    #Df[i,4]=((m-1)/(m))*Df[i,4]
    Df2=copy(df2)
    A32,A42=XiMass(Df,Df2);
    V=pinv(A32)*A42-pinv(A31)*A41
    
    for j in 1:m
        #=
        Df=copy(df)
        Df[i,4]=((m-j)/(m))*Df[i,4]
        #Df[i,4]=((m-1)/(m))*Df[i,4]
        Df2=copy(df2)
        A32,A42=XiMass(Df,Df2);
        V=pinv(A32)*A42-pinv(A31)*A41
        =#
        λper=λ+diag(IQ*V*Q)
        Δper = Diagonal(λper)
        A=reshape(repeat(transpose.(λ), outer = 2*n),(2*n,2*n))
        AT=transpose(A)
        DAT=inv.(AT-A+I)-I
        D=(IQ*V*Q).*DAT
        Qpert=Q+Q*D
        #A31,A41=A32,A42
        λ=λper
        Q=Qpert
        IQ=pinv(Qpert,0.03)
    end
    Δ = Diagonal(λ)
    A31,A41=XiMass(df,df2);
    #f1 = (t) -> Q*exp(t*Δ)*pinv(Q)*Y0 + Q*((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-1],1)))+Diagonal(vcat(zeros(2*n-1),t)))*pinv(Q)*(pinv(A31)*b)
    #f1 = (t) -> Q*exp(t*Δ)*IQ*Y0 + Q*((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-1],1)))+Diagonal(vcat(zeros(2*n-1),t)))*IQ*(pinv(A31)*b)
    nozev=sum(abs.(λ).<1e-10)
    f1 = (t) -> Q*exp(t*Δ)*IQ*Y0 + Q*((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-nozev],fill(1,nozev))))+Diagonal(vcat(zeros(2*n-nozev),fill(t,nozev))))*IQ*(pinv(A31)*b)
    

    #=
    #Single Phase
    Df=copy(df)
    Df[i,4]=(2/3)*Df[i,4]
    Df2=copy(df2)
    A32,A42=XiMass(Df,Df2);
    V=pinv(A32)*A42-pinv(A31)*A41
    λper=λ+diag(IQ*V*Q)
    Δper = Diagonal(λper)

    A=reshape(repeat(transpose.(λ), outer = 2*n),(2*n,2*n))
    AT=transpose(A)
    DAT=inv.(AT-A+I)-I
    D=(IQ*V*Q).*DAT
    Qpert=Q+Q*D

    #Two Phase
    Cf=copy(df)
    Cf[i,4]=(1/3)*Cf[i,4]
    Cf2=copy(df2)
    A34,A44=XiMass(Cf,Cf2);
    V=pinv(A34)*A44-pinv(A32)*A42
    λper=λper+diag(pinv(Qpert)*V*Qpert)
    Δper = Diagonal(λper)

    A=reshape(repeat(transpose.(λ), outer = 2*n),(2*n,2*n))
    AT=transpose(A)
    DAT=inv.(AT-A+I)-I
    D=(pinv(Qpert)*V*Qpert).*DAT
    Qpert=Qpert+Qpert*D
    

    #Three Phase
    Bf=copy(df)
    Bf[i,4]=0*Bf[i,4]
    Bf2=copy(df2)
    A33,A43=XiMass(Bf,Bf2);
    V1=pinv(A33)*A43-pinv(A32)*A42
    λper2=λper+diag(pinv(Qpert)*V1*Qpert)
    Δper2 = Diagonal(λper2)

    
    A2=reshape(repeat(transpose.(λper), outer = 2*n),(2*n,2*n))
    AT2=transpose(A2)
    DAT2=inv.(AT2-A2+I)-I
    D=(pinv(Qpert)*V1*Qpert).*DAT2
    Qpert2=Qpert+Qpert*D
    
    f1 = (t) -> Qpert2*exp(t*Δper2)*pinv(Qpert2)*Y0 + Qpert2*((exp(t*Δper2)-I).*pinv(Diagonal(vcat(λper2[1:end-1],1)))+Diagonal(vcat(zeros(2*n-1),t)))*pinv(Qpert2)*(pinv(A31)*b)
    =#
    
    return f1
end


function Ysol000Nunif(df,df2,i,m)
    n=df.Nnodes[1]; #Number of of buses    
    P1=df2.PowerInjections; #Vector of power injections 1
    L1=Laplatian(df); #Laplacian Matrix F
    
    Theta0=pinv(L1)*P1; #Steady-state solution
    Y0=[zeros(n);Theta0];
    b=[P1;zeros(n)];

    
    #NoFault
    A31,A41=XiMass(df,df2);
    XIEig1=eigen(pinv(A31)*A41);
    Q1 = XIEig1.vectors;
    λ  = XIEig1.values;
    Δ = Diagonal(λ);
    
    Q=Q1
    #IQ=pinv(Q)
    IQ=pinv(Q)
    
    #h=reverse(collect(range(0, stop = 2/3, length = m)))
    #p=((m-1)/(m))
    #h=pushfirst!(h,p)

    w=0.79
    x=1-w
    h=reverse(collect(range(0, stop = 2/3, length = trunc(Int, m*w))))
    p=((m-1)/(m))
    h2=reverse(collect(range(2/3, stop=p, length = trunc(Int, m*x))))
    h=unique!(vcat(h2,h))

    #=
    Df=copy(df)
    Df[i,4]=((m-1)/(m))*Df[i,4]
    Df2=copy(df2)
    A32,A42=XiMass(Df,Df2);
    V=pinv(A32)*A42-pinv(A31)*A41
    =#

    for j in h

        
        Df=copy(df)
        Df[i,4]=j*Df[i,4]
        Df2=copy(df2)
        A32,A42=XiMass(Df,Df2);
        V=pinv(A32)*A42-pinv(A31)*A41

        
        λper=λ+diag(IQ*V*Q)
        #IQ=W0,V=A1,Q=V0
        #λ+diag(W0*A1*V0)
        Δper = Diagonal(λper)
        A=reshape(repeat(transpose.(λ), outer = 2*n),(2*n,2*n))
        AT=transpose(A)
        DAT=inv.(AT-A+I)-I
        D=DAT.*(IQ*V*Q)
        Qpert=Q+Q*D
        #
        #
        #
        #
        A31,A41=A32,A42
        λ=λper
        Q=Qpert
        #IQ=pinv(Qpert)
        IQ=pinv(Qpert,0.03)
    end
    Δ = Diagonal(λ)
    A31,A41=XiMass(df,df2);
    nozev=sum(abs.(λ).<1e-10)
    f1 = (t) -> Q*exp(t*Δ)*IQ*Y0 + Q*((exp(t*Δ)-I).*pinv(Diagonal(vcat(λ[1:end-nozev],fill(1,nozev))))+Diagonal(vcat(zeros(2*n-nozev),fill(t,nozev))))*IQ*(pinv(A31)*b)

    return f1
end



function AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,value=0.006)
    #AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,value=0.006)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    f1=Ysol0(df,df2,Df,Df2)
    Tempo1=[T1:value:T2;]
    Tempo2=[T2:value:T3;]
    RES=hcat(real.(f1.(Tempo1))...)
    #for i in 2:length(Tempo1)
    #   RES=hcat(RES,f1(Tempo1[i])) 
    #end

    
    Y1=real(f1(T2))
    f2=Ysol2(df,df2,Df,Df2,Y1)
    #for i in 2:length(Tempo2)
    #   RES=hcat(RES,f2(Tempo2[i])) 
    #end
    RES=hcat(RES,hcat(real.(f2.(Tempo2.-T2))...))
    MAT=RES
    J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    DFX.Time = range(0, T3, length=size(DFX)[1])
    #DFdX=DataFrame(J[:,1:n],:auto)
    #DFdX.Time = range(0, T3, length=size(DFdX)[1])
    return DFX
end

function AnalyticalSolutionApprox0(df,df2,Df,Df2,Bf,Bf2,T1,T2,T3,value=0.006)
    #AnalyticalSolutionApprox0(df,df2,Df,Df2,Bf,Bf2,T1,T2,T3,value=0.006)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    f1=Ysol00(df,df2,Df,Df2,Bf,Bf2)
    Tempo1=[T1:value:T2;]
    Tempo2=[T2:value:T3;]
    RES=hcat(real.(f1.(Tempo1))...)
    #for i in 2:length(Tempo1)
    #   RES=hcat(RES,f1(Tempo1[i])) 
    #end

    
    Y1=real(f1(T2))
    f2=Ysol2(df,df2,Df,Df2,Y1)
    #for i in 2:length(Tempo2)
    #   RES=hcat(RES,f2(Tempo2[i])) 
    #end
    RES=hcat(RES,hcat(real.(f2.(Tempo2.-T2))...))
    MAT=RES
    J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    DFX.Time = range(0, T3, length=size(DFX)[1])
    #DFdX=DataFrame(J[:,1:n],:auto)
    #DFdX.Time = range(0, T3, length=size(DFdX)[1])
    return DFX
end

function AnalyticalSolutionApprox00(df,df2,Df,Df2,T1,T2,T3,i,m,value=0.006)
    # AnalyticalSolutionApprox00(df,df2,Df,Df2,T1,T2,T3,i,m,value=0.006)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    f1=Ysol000(df,df2,i,m)
    Tempo1=[T1:value:T2;]
    Tempo2=[T2:value:T3;]
    RES=hcat(real.(f1.(Tempo1))...)
    #for i in 2:length(Tempo1)
    #   RES=hcat(RES,f1(Tempo1[i])) 
    #end

    
    Y1=real(f1(T2))
    f2=Ysol2(df,df2,Df,Df2,Y1)
    #for i in 2:length(Tempo2)
    #   RES=hcat(RES,f2(Tempo2[i])) 
    #end
    RES=hcat(RES,hcat(real.(f2.(Tempo2.-T2))...))
    MAT=RES
    J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    DFX.Time = range(0, T3, length=size(DFX)[1])
    #DFdX=DataFrame(J[:,1:n],:auto)
    #DFdX.Time = range(0, T3, length=size(DFdX)[1])
    return DFX
end

function AnalyticalSolutionApprox00Nunif(df,df2,Df,Df2,T1,T2,T3,i,m,value=0.006)
    f1=Ysol000Nunif(df,df2,i,m)
    Tempo1=[T1:value:T2;]
    Tempo2=[T2:value:T3;]
    RES=hcat(real.(f1.(Tempo1))...)
    
    Y1=real(f1(T2))
    f2=Ysol2(df,df2,Df,Df2,Y1)
    
    RES=hcat(RES,hcat(real.(f2.(Tempo2.-T2))...))
    MAT=RES
    J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    DFX.Time = range(0, T3, length=size(DFX)[1])
    return DFX
end


function AnalyticalSolutionDX(df,df2,Df,Df2,T1,T2,T3,value=0.006)
    #AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,value=0.003)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    f1,Qpert,Δpert=Ysol1(df,df2,Df,Df2)
    Tempo1=[T1:value:T2;]
    Tempo2=[T2:value:T3;]
    RES=hcat(real.(f1.(Tempo1))...)
    #for i in 2:length(Tempo1)
    #   RES=hcat(RES,f1(Tempo1[i])) 
    #end
    Y1=real(f1(T2))
    f2=Ysol2(df,df2,Df,Df2,Y1)
    #for i in 2:length(Tempo2)
    #   RES=hcat(RES,f2(Tempo2[i])) 
    #end
    RES=hcat(RES,hcat(real.(f2.(Tempo2.-T2))...))
    MAT=RES
    J=transpose(MAT)
    #DFX=DataFrame(J[:,n+1:2*n],:auto)
    #DFX.Time = range(0, T3, length=size(DFX)[1])
    DFdX=DataFrame(J[:,1:n],:auto)
    DFdX.Time = range(0, T3, length=size(DFdX)[1])
    return DFdX
end

function AnalyticalSolutionApproxDX(df,df2,Df,Df2,T1,T2,T3,value=0.006)
    #AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,value=0.01)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    f1=Ysol0(df,df2,Df,Df2)
    Tempo1=[T1:value:T2;]
    Tempo2=[T2:value:T3;]
    RES=hcat(real.(f1.(Tempo1))...)
    #for i in 2:length(Tempo1)
    #   RES=hcat(RES,f1(Tempo1[i])) 
    #end

    
    Y1=real(f1(T2))
    f2=Ysol2(df,df2,Df,Df2,Y1)
    #for i in 2:length(Tempo2)
    #   RES=hcat(RES,f2(Tempo2[i])) 
    #end
    RES=hcat(RES,hcat(real.(f2.(Tempo2.-T2))...))
    MAT=RES
    J=transpose(MAT)
    #DFX=DataFrame(J[:,n+1:2*n],:auto)
    #DFX.Time = range(0, T3, length=size(DFX)[1])
    DFdX=DataFrame(J[:,1:n],:auto)
    DFdX.Time = range(0, T3, length=size(DFdX)[1])
    return DFdX
end


function SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi=0.0,value=0.001)
#function SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi=0.01,value=0.002)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    n=df.Nnodes[1] #Number of of buses    
    P1=df2.PowerInjections #Vector of power injections 1
    A3,A4=XiMass(df,df2) #Xi and Mass Matrix 1
    P2=Df2.PowerInjections #Vector of power injections 2s
    A5,A6=XiMass(Df,Df2) #Xi and Mass Matrix 2    
    function f(du,u,p,t)
        if (t>=T1)&(t<=T2)
            du.=A6*u+[P2;zeros(n)]
            #du.=A5*u+[P2;zeros(n)]
        else
            du.=A4*u+[P1;zeros(n)]
            #du.=A3*u+[P1;zeros(n)]
        end
    end
    function noise(du, u, p, t)
        for i in 1:n
            du[i] = noi
        end 
    end
    L1=Laplatian(df) # Laplatian matrix 1
    
    Theta0=pinv(L1)*P1 #Steady-state solution
    
    #Theta0=df2.Theta0

    u0 = [zeros(n);Theta0];
    
    tspan1 = (0.0,T3);
    ###
    #SDE
    prob1 = SDEProblem(SDEFunction(f,noise; mass_matrix=A5),u0,tspan1,noise);
    #prob1 = SDEProblem(SDEFunction(f,noise; mass_matrix=A4),u0,tspan1,noise);
    sol1 = solve(prob1,ImplicitEM(),dtmax=value);
    
    #sol1 = solve(prob1,ImplicitEulerHeun(),dtmax=value);
    
    #sol1 = solve(prob1,EM(),dt=value,dtmax=value);
    #sol1 = solve(prob1,SRA3(),dtmax=value);
    #RES = sol1.u[1];
    #for i in 1:length(sol1.t)
    #   RES=hcat(RES,sol1.u[i]) 
    #end
    #RES=hcat(RES,hcat(sol1.u...))
    RES=hcat(sol1.u...)
    MAT=RES
    J=transpose(MAT)    
    
    #Store Info as DataFrame
    #MAT=RES[:,2:end]
    #J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    #DFdX=DataFrame(J[:,1:n],:auto)
    #=
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:n#size(DFX)[2]
        for j in 1:n#size(DFX)[2]
            if j>i
                I=string(i)
                J=string(j)
                comma=","
                str="X"*I*comma*J
                df3=DataFrame(str=>DFX[:,i]-DFX[:,j])
                DFX=hcat(DFX,df3)
                                                
            end
        end
    end
    =#
    DFX.Time = range(0, T3, length=size(DFX)[1])
    #CSV.write("DFX1.csv", DFX)
    #DFdX.Time = range(0, T3, length=size(DFdX)[1])
    #CSV.write("DFdX.csv", DFdX)
    return DFX
end

function SimulationDEFSDE3(df,df2,Df,Df2,T1,T2,T3,noi=0.01,value=0.001)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    n=df.Nnodes[1] #Number of of buses    
    P1=df2.PowerInjections #Vector of power injections 1
    A3,A4=XiMass(df,df2) #Xi and Mass Matrix 1
    P2=Df2.PowerInjections #Vector of power injections 2s
    A5,A6=XiMass(Df,Df2) #Xi and Mass Matrix 2    
    function f(du,u,p,t)
        if (t>=T1)&(t<=T2)
            du.=A6*u+[P2;zeros(n)]
            #du.=A5*u+[P2;zeros(n)]
        else
            du.=A4*u+[P1;zeros(n)]
            #du.=A3*u+[P1;zeros(n)]
        end
    end
    function noise(du, u, p, t)
        for i in 1:n
            du[i] = noi
        end 
    end
    L1=Laplatian(df) # Laplatian matrix 1
    
    Theta0=pinv(L1)*P1 #Steady-state solution
    
    #Theta0=df2.Theta0

    u0 = [zeros(n);Theta0];
    
    tspan1 = (0.0,T3);
    ###
    #SDE
    prob1 = SDEProblem(SDEFunction(f,noise; mass_matrix=A5),u0,tspan1,noise);
    #prob1 = SDEProblem(SDEFunction(f,noise; mass_matrix=A4),u0,tspan1,noise);
    #sol1 = solve(prob1,ImplicitEM(),dtmax=value);
    sol1 = solve(prob1,ImplicitEulerHeun(),dtmax=value);
    #sol1 = solve(prob1,EM(),dt=value,dtmax=value);
    #sol1 = solve(prob1,SRA3(),dtmax=value);
    
    #RES = sol1.u[1];
    #for i in 1:length(sol1.t)
    #   RES=hcat(RES,sol1.u[i]) 
    #end

    #RES=hcat(RES,hcat(sol1.u...))
    #MAT=RES[:,2:end]
    RES=hcat(sol1.u...)
    MAT=RES
    J=transpose(MAT)    
    
    #Store Info as DataFrame
    #MAT=RES[:,2:end]
    #J=transpose(MAT)
    #DFX=DataFrame(J[:,n+1:2*n],:auto)
    DFdX=DataFrame(J[:,1:n],:auto)
    #=
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:n#size(DFX)[2]
        for j in 1:n#size(DFX)[2]
            if j>i
                I=string(i)
                J=string(j)
                comma=","
                str="X"*I*comma*J
                df3=DataFrame(str=>DFX[:,i]-DFX[:,j])
                DFX=hcat(DFX,df3)
                                                
            end
        end
    end
    =#
    #DFX.Time = range(0, T3, length=size(DFX)[1])
    #CSV.write("DFX1.csv", DFX)
    DFdX.Time = range(0, T3, length=size(DFdX)[1])
    #CSV.write("DFdX.csv", DFdX)
    return DFdX
end

function FreqAnalysis(df,df2,Df,Df2,T1,T3,TN=2.0,noi=0.01,value=0.002,lim=1.5)
    T=0.5:0.1:TN
    Fmax=[0.0]
    h=0
    for t in T
        T2=T1+t
        H=SimulationDEFSDE3(df,df2,Df,Df2,T1,T2,T3,0.0)
        p=stack(H[:,1:(end-1)])
        f=maximum(abs.(p[:,2]./(2*pi)))
        Fmax=push!(Fmax,f)
        h=t
        if any(Fmax.>=lim)
            break
        else 
            continue
        end 
    end 
    A=h
    B=Fmax[end]
    return A,B
end






function Analysis0(DFX,df,df2)
    th=df.ThetaMax
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str="X"*I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+2:size(DFX)[2]]
    B=Bool(1)
    V=Bool(0)
    mst=zeros(size(DF)[1], 1)
    for i in 1:N
        h=Bool.(abs.(DF[:,i]).<=th[i])
        mst=hcat(mst,h)
        B=B.&&h
    end
    #Outside
    tra=findall(mst[:,Not([1])].==0)
    TRAC=hcat(getindex.(tra, 1), getindex.(tra,2))
    TRAC=TRAC[sortperm(TRAC[:, 1]), :]
    if size(TRAC)[1]>0#sum(V)>0
        # First Exit Time
        Time1 = DFX.Time[minimum(TRAC[:,1])]
        #Time1 = minimum(DFX.Time[V])
        alpha = DFX.Time.>Time1
        beta = alpha.&&B
        if sum(beta)>0
            # First Return Time
            Time2 = minimum(DFX.Time[beta])
        else 
            Time2 = "Undefined"
        end
        if sum(Bool.(pushfirst!(diff(TRAC[:,1]),2).>1))>0
            ExitTimes=DFX.Time[TRAC[:, 1][Bool.(pushfirst!(diff(TRAC[:,1]),2).>1)]]
        else
            ExitTimes=[Time1]
        end
        g=(DFX.Time .∉ Ref(ExitTimes)).*(1:size(DFX.Time)[1])
        if sum(Bool.(pushfirst!(diff(g),2).>1))>0
            EntryTimes=DFX.Time[Bool.(pushfirst!(diff(g),2).>1)]
            # Last Return Time
            Time3=EntryTimes[end]
        else
            EntryTimes=[Time2]
            Time3=Time2
        end
    else
        Time1 = NaN
        Time2 = NaN
        Time3 = NaN
        ExitTimes = [NaN]
        EntryTimes = [NaN]
    end
    Results = DataFrame(FirstExitTime=Float64[],
        FirstReturnTime=Float64[],
        LastReturnTime=Float64[])
    push!(Results,(Time1,Time2,Time3)) 
    return Results
end



function Analysis(DFX,df,df2)
    th=df.ThetaMax
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str="X"*I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+2:size(DFX)[2]]
    B=Bool(1)
    V=Bool(0)
    #vic=[]
    mst=zeros(size(DF)[1], 1)
    for i in 1:N
        h=Bool.(abs.(DF[:,i]).<=th[i])
        #push!(vic,!prod(Bool.(h)))
        mst=hcat(mst,h)
        B=B.&&h
        #j=Bool.(abs.(DF[:,i]).>th[i])
        #V=V.||j
    end
    #Outside
    tra=findall(mst[:,Not([1])].==0)
    TRAC=hcat(getindex.(tra, 1), getindex.(tra,2))
    TRAC=TRAC[sortperm(TRAC[:, 1]), :]
    #vic=Bool.(vic)
    #Mean1=mean(B)
    #Mean2 = tspan[2]*Mean1
    #Mean3 = tspan[2]-Mean2
    if size(TRAC)[1]>0#sum(V)>0
        # First Exit Time
        Time1 = DFX.Time[minimum(TRAC[:,1])]
        #Time1 = minimum(DFX.Time[V])
        alpha = DFX.Time.>Time1
        beta = alpha.&&B
        if sum(beta)>0
            # First Return Time
            Time2 = minimum(DFX.Time[beta])
        else 
            Time2 = "Undefined"
        end
        if sum(Bool.(pushfirst!(diff(TRAC[:,1]),2).>1))>0
            ExitTimes=DFX.Time[TRAC[:, 1][Bool.(pushfirst!(diff(TRAC[:,1]),2).>1)]]
        else
            ExitTimes=[Time1]
        end
        g=(DFX.Time .∉ Ref(ExitTimes)).*(1:size(DFX.Time)[1])
        if sum(Bool.(pushfirst!(diff(g),2).>1))>0
            EntryTimes=DFX.Time[Bool.(pushfirst!(diff(g),2).>1)]
        else
            EntryTimes=[Time2]
        end
        
        TRAC2=unique(TRAC[:,2])
        #FL=join(names(DF)[TRAC2],"|")
        FL2=names(DF)[TRAC2]
        FFL=names(DF)[TRAC2[1]]
        #if size(TRAC)[1]>0#sum(vic)>0
        #    TRAC2=unique(TRAC[:,2])
        #    FL=join(names(DF)[TRAC2],"|")
        #else
        #    FL = "Undefined"
        #    FFL = "Undefined"
        #end 
    else
        Time1 = NaN
        Time2 = NaN
        #FL = "Undefined"
        FFL = "Undefined"
        ExitTimes = [NaN]
        EntryTimes = [NaN]
        FL2 = ["Undefined"]
        #OutTimes = []
        #InTimes = []
    end
    Results = DataFrame(FirstExitTime=Float64[],
        #Union{Nothing, Float64}[]
        #
        FirstReturnTime=Float64[],
        #FirstReturnTime=Any[],
        FirstFaultLine=String[],
        #FaultyLines=Array{String,1}[],
        ExitTimes=Array{Float64,1}[],
        EntryTimes=Array{Float64,1}[],
        FaultyLines2=Array{String,1}[])
        #MeanIns = Mean2,
        #MeanOut = Mean3, 
    push!(Results,(Time1,Time2,FFL,ExitTimes,EntryTimes,FL2)) 
    return Results
end


function AnalysisBeta(DFX,df,df2)
    th=df.ThetaMax
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str="X"*I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+1:size(DFX)[2]]
    TDF=DataFrame([[names(DF)]; collect.(eachrow(DF))], [:column; Symbol.(axes(DF, 1))])
    ExitTimes=[]
    EntryTimes=[]
    FL2=[]
    PhaseDifferencesTransposed=TDF[2:end,2:end]
    for i in 2:(ncol(TDF[2:end,2:end]))
        is_inside0 = !any(abs.(PhaseDifferencesTransposed[:,i-1]).>th)
        is_inside1 = !any(abs.(PhaseDifferencesTransposed[:,i]).>th)
        if is_inside0 && !is_inside1
            push!(ExitTimes,DF.Time[i])
            push!(FL2,TDF.column[2:end][(abs.(PhaseDifferencesTransposed[:,i]).>th)])
        end
        if is_inside1 && !is_inside0 
            push!(EntryTimes,DF.Time[i])
        end
        
    end
    if size(ExitTimes)[1]>0
        Time1 = ExitTimes[1]
    else
        Time1 = NaN
        ExitTimes = [NaN]
    end
    if size(EntryTimes)[1]>0
        Time2 = EntryTimes[1]
    else
        Time2 = NaN
        EntryTimes = [NaN]
    end
    if size(FL2)[1]>0
        FFL = FL2[1][1]
    else
        FFL = "Undefined"
        FL2 = ["Undefined"]
    end
    Results = DataFrame(FirstExitTime=Float64[],
        
        FirstReturnTime=Float64[],
        
        FirstFaultLines=Any[],
        
        ExitTimes=Array{Float64,1}[],
        EntryTimes=Array{Float64,1}[])
        #,AllFaultyLinesEx=Any[])
        
    push!(Results,(Time1,Time2,FFL,
                ExitTimes,EntryTimes))
            #FL2)) 
    return Results
end

function AnalysisDelta(DFX,df,df2)
    th=df.ThetaMax
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str="X"*I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+1:size(DFX)[2]]
    TDF=DataFrame([[names(DF)]; collect.(eachrow(DF))], [:column; Symbol.(axes(DF, 1))])    
    tdf=copy(TDF)
    PhaseDifferencesTransposed=TDF[2:end,2:end]
    OutsideTimes=[]
    FL2=[]
    for i in 2:(ncol(TDF[2:end,2:end]))
        is_outside = any(abs.(PhaseDifferencesTransposed[:,i]).>th)
        if is_outside
            push!(OutsideTimes,DF.Time[i])
            push!(FL2,TDF.column[2:end][(abs.(PhaseDifferencesTransposed[:,i]).>th)])
        end    
    end
    if size(OutsideTimes)[1]>0
        OutsideTimes = OutsideTimes
    else
        OutsideTimes = [NaN]
    end
    Results = DataFrame(OutsideTimes=Any[],
        AllFaultyLinesEx=Any[])
        
    push!(Results,(OutsideTimes,FL2)) 
    return Results    
end

function AnalysisGamma(DFX,df,df2)
    th=df.ThetaMax
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str="X"*I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+1:size(DFX)[2]]
    TDF=DataFrame([[names(DF)]; collect.(eachrow(DF))], [:column; Symbol.(axes(DF, 1))])    
    tdf=copy(TDF)
    #PhaseDifferencesTransposed=TDF[2:end,2:end]
    #for i in 1:(ncol(TDF[2:end,2:end]))
    #    tdf[2:end,i+1]=(abs.(PhaseDifferencesTransposed[:,i]).>th)    
    #end    
    tdf[2:end,2:end]=(abs.(TDF[2:end,2:end]).>th)
    return tdf
end

function AnalysisGammaIndv(DFX,df,df2)
    th=df.ThetaMax
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str=I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+1:size(DFX)[2]]
    TDF=DataFrame([[names(DF)]; collect.(eachrow(DF))], [:column; Symbol.(axes(DF, 1))])    
    tdf=copy(TDF)
    #PhaseDifferencesTransposed=TDF[2:end,2:end]
    #for i in 1:(ncol(TDF[2:end,2:end]))
    #    tdf[2:end,i+1]=(abs.(PhaseDifferencesTransposed[:,i]).>th)    
    #end    
    tdf[2:end,2:end]=(abs.(TDF[2:end,2:end]).>th)
    return tdf
end

function AnalysisGamma2(DFX,df,df2)
    th=df.ThetaMax
    B=df.Susceptance
    F=df.F
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str="X"*I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+1:size(DFX)[2]]
    TDF=DataFrame([[names(DF)]; collect.(eachrow(DF))], [:column; Symbol.(axes(DF, 1))])    
    tdf=copy(TDF)
    #display(tdf)
    
    #PhaseDifferencesTransposed=TDF[2:end,2:end]
    
    #for i in 1:(ncol(TDF[2:end,2:end]))
    #    tdf[2:end,i+1]=(abs.(PhaseDifferencesTransposed[:,i])./th)    
    #    #tdf[2:end,i+1]=(B.*abs.(PhaseDifferencesTransposed[:,i])./F)    
    #end    
    tdf[2:end,2:end]=(abs.(TDF[2:end,2:end])./th)

    return tdf
end


function AnalysisGamma3(DFX,df,df2)
    th=df.ThetaMax
    B=df.Susceptance
    F=df.F
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str="X"*I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+1:size(DFX)[2]]
    TDF=DataFrame([[names(DF)]; collect.(eachrow(DF))], [:column; Symbol.(axes(DF, 1))])    
    tdf=copy(TDF)
    #display(tdf)
    
    PhaseDifferencesTransposed=TDF[2:end,2:end]
    
    #for i in 1:(ncol(TDF[2:end,2:end]))
    #    tdf[2:end,i+1]=PhaseDifferencesTransposed[:,i]./th    
    #    #tdf[2:end,i+1]=(B.*abs.(PhaseDifferencesTransposed[:,i])./F)        
    #end    
    tdf[2:end,2:end]=(TDF[2:end,2:end].>th)

    return tdf
end


#=
function OverheatingIndicator(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2)
    Opt=tdf[2:end,2:end]
    delta=values(tdf[1,2:end])[2]-values(tdf[1,2:end])[1]
    j=Opt.*delta
    j_selected=j[:,sum.(eachcol(j)).>0]
    arr=filter(x -> x > 0, vcat(eachcol(j_selected)...))
    if size(arr)[1]==0
        return(0.0)
    else
        return(sum(arr))
    end 
end

=#

#=
function OverheatingIndicator(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2);
    j_selected=(tdf[2:end,2:end].*(tdf[1,3]-tdf[1,2]))[:,sum.(eachcol(tdf[2:end,2:end].*(tdf[1,3]-tdf[1,2]))).>0];
    arr=filter(x -> x > 0, vcat(eachcol(j_selected)...));
    if size(arr)[1]==0
        return(0.0)
    else
        return(sum(arr))
    end 
end
=#


#=
function OverheatingIndicator(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2);
    M = Matrix(tdf[2:end, 2:end])
    Δ = tdf[1,3] - tdf[1,2]

    # Compute column sums without building M .* Δ
    colsum = vec(sum(M, dims=1)) .* Δ
    keep = colsum .> 0

    # No column kept?
    any(keep) || return 0.0

    # Sum positives directly
    s = 0.0
    for j in eachindex(keep)
        keep[j] || continue
        @inbounds for x in @view(M[:,j])
            v = x * Δ
            v > 0 && (s += v)
        end
    end
    return s 
end
=#

function OverheatingIndicator(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2);
    M = Float32.(Matrix(tdf[2:end, 2:end]))
    Δ = Float32(tdf[1,3] - tdf[1,2])

    colsum = vec(sum(M, dims=1)) .* Δ
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




function OverheatingIndicatorIndv(DFX,df,df2)
    tdf=AnalysisGammaIndv(DFX,df,df2)
    delta=values(tdf[1,2:end])[2]-values(tdf[1,2:end])[1]
    Opt=tdf[2:end,1:end]
    Opt[1:end,2:end]=Opt[1:end,2:end].*delta
    Opt.S = sum.(eachrow(Opt[:, names(Opt, Real)]))
    comma=","
    sim="i"*comma*"j"
    rename!(Opt,:column => sim)
    rename!(Opt,:S => :Sij)
    disc=Opt[:,[1,end]]
    return(disc)
end


function OverheatingIndicatorREDEF(DFX,df,df2)
    tdf=AnalysisGamma2(DFX,df,df2)
    delta=values(tdf[1,2:end])[2]-values(tdf[1,2:end])[1]
    tdf[2:end,2:end]=Float64.(tdf[2:end,2:end].>=1).*tdf[2:end,2:end]
    h=tdf[2:end,2:end]
    h_selected=h[:,sum.(eachcol(h)).>0]
    arr=filter(x -> x > 0, vcat(eachcol(h_selected)...))
    if size(arr)[1]==0
        return(0.0)
    else 
        arrpro=1 ./((4 ./ (arr .- 1)).^(2.86))
        return(sum(delta.*arrpro))
    end 
end

function OverheatingIndicator2(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2)
    Opt=tdf[2:end,2:end]
    delta=values(tdf[1,2:end])[2]-values(tdf[1,2:end])[1]
    Arr=[]
    for i in 1:(size(Opt)[2])
        push!(Arr,delta*sum(Opt[:,i]))
    end
    str1="Time"
    str2="Stp"
    df1=DataFrame(str1=>collect(tdf[1,2:end]))
    df2=DataFrame(str2=>cumsum(Arr))
    Vic=hcat(df1,df2)
    return(Vic)
end


#=
function AnalysisDelta(DFX,df,df2)
    th=df.ThetaMax
    n=df.Nnodes[1]
    N=nrow(df)
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:N
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        str="X"*I*comma*J
        df3=DataFrame(str=>DFX[:,df.From[i]]-DFX[:,df.To[i]])
        DFX=hcat(DFX,df3)
    end
    DF=DFX[:,n+1:size(DFX)[2]]
    TDF=DataFrame([[names(DF)]; collect.(eachrow(DF))], [:column; Symbol.(axes(DF, 1))])    
    tdf=copy(TDF)
    PhaseDifferencesTransposed=TDF[2:end,2:end]
    OutsideTimes=[]
    FL2=[]
    for i in 2:(ncol(TDF[2:end,2:end]))
        is_outside = any(abs.(PhaseDifferencesTransposed[:,i]).>th)
        if is_outside
            push!(OutsideTimes,DF.Time[i])
            push!(FL2,TDF.column[2:end][(abs.(PhaseDifferencesTransposed[:,i]).>th)])
        end    
    end
    if size(OutsideTimes)[1]>0
        OutsideTimes = OutsideTimes
    else
        OutsideTimes = [NaN]
    end
    Results = DataFrame(OutsideTimes=Any[],
        AllFaultyLinesEx=Any[])
        
    push!(Results,(OutsideTimes,FL2)) 
    return Results    
end
=#

function AnalysisEpsilon(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2)[2:end,2:end]
    Time=collect(values(AnalysisGamma(DFX,df,df2)[1,2:end]))
    Lines=AnalysisGamma(DFX,df,df2)[2:end,1]
    Results = DataFrame(TimesOfChange=Any[],
    FaultyLines=Any[])
    for i in 1:(ncol(tdf)-1)
        D=any(tdf[1:end,i+1].!=tdf[1:end,i])
        if D
            push!(Results,(Time[i+1],Lines[Bool.(tdf[1:end,i+1])]))
        end
    end
    return(Results)
end

function AnalysisZeta(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2)
    G=permutedims(tdf, 1)[:,3:end]
    Time=collect(values(AnalysisGamma(DFX,df,df2)[1,2:end]))
    Lines=AnalysisGamma(DFX,df,df2)[2:end,1]
    Results = DataFrame(FaultyLines=Any[],
                        TimesOfChange=Any[])
    for i in 1:(ncol(G))
        D=any(G[1:end,i].==1.0)
        if D
            B=permutedims(tdf, 1)[1:end,2:end][Bool.(permutedims(tdf, 1)[1:end,i+2].==1.0),:]
            Time=[values(B[:,1])[1]]
            #Time=[]
            for j in 1:(nrow(B)-1)
                E=any(values(B[j+1,2:end]).!=values(B[j,2:end]))
                if E
                    push!(Time,values(B[:,1])[j+1])
                end
            end
            push!(Results,(Lines[i],Time))
        end
    end
    return(Results)
end

function AnalysisEta(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2)
    G=permutedims(tdf, 1)[:,3:end]
    Tim=collect(values(AnalysisGamma(DFX,df,df2)[1,2:end]))
    Lines=AnalysisGamma(DFX,df,df2)[2:end,1]
    Results = DataFrame(FaultyLines=Any[],
                        TimesOfChange=Any[])
    for i in 1:(ncol(G))
        h=findall(G[:,i].==1.0)
        if size(h)[1]>=1
            Time=[Tim[h[1]]]
            for j in 1:(size(h)[1]-1)
                if abs(h[j+1]-h[j])>1
                    push!(Time,Tim[h[j]+1])
                    push!(Time,Tim[h[j+1]])
                end
            end
            push!(Time,Tim[h[end]+1])
            push!(Results,(Lines[i],Time))
        end
    end
    return(Results)
end

function AnalysisTheta(DFX,df,df2)
    tdf=AnalysisGamma(DFX,df,df2)
    G=permutedims(tdf, 1)[:,3:end]
    Tim=collect(values(AnalysisGamma(DFX,df,df2)[1,2:end]))
    Lines=AnalysisGamma(DFX,df,df2)[2:end,1]
    Results = DataFrame(FaultyLines=Any[],
                        TimesOfChange=Any[])
    for i in 1:(ncol(G))
        h=findall(G[:,i].==1.0)
        if size(h)[1]>=1
            Time=[Tim[h[1]]]
            for j in 1:(size(h)[1]-1)
                if abs(h[j+1]-h[j])>1
                    push!(Time,Tim[h[j]+1])
                    push!(Time,Tim[h[j+1]])
                end
            end
            push!(Time,Tim[h[end]+1])
            push!(Results,(Lines[i],Time))
        end
    end
    return(Results)
end

function ComparisonPhase(df,df2,T1,T2,T3,alpha)
    G=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)
        Astacked_df = [stack(ADFX[:,1:(end-1)])][1]
        Time=repeat(ADFX[:,end],Int(size(Astacked_df)[1]/size(ADFX)[1]))
        Astacked_df.Time = Time
        Analytical=Astacked_df
    
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)
        Bstacked_df = [stack(BDFX[:,1:(end-1)])][1]
        Time=repeat(BDFX[:,end],Int(size(Bstacked_df)[1]/size(BDFX)[1]))
        Bstacked_df.Time = Time
        Approx=Bstacked_df
        #push!(G,maximum(abs.(Approx[:,2]-Analytical[:,2])))
        push!(G, sqrt((Approx[:,3][2])*sum((Approx[:,2]-Analytical[:,2]).^2)) )
    end 
    
    return maximum(G)
    
end 

function ComparisonFrequency(df,df2,T1,T2,T3,alpha)
    G=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)
        Astacked_df = [stack(ADFdX[:,1:(end-1)])][1]
        Time=repeat(ADFdX[:,end],Int(size(Astacked_df)[1]/size(ADFdX)[1]))
        Astacked_df.Time = Time
        Analytical=Astacked_df
    
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)
        Bstacked_df = [stack(BDFdX[:,1:(end-1)])][1]
        Time=repeat(BDFdX[:,end],Int(size(Bstacked_df)[1]/size(BDFdX)[1]))
        Bstacked_df.Time = Time
        Approx=Bstacked_df
        #push!(G,maximum(abs.(Approx[:,2]-Analytical[:,2])))
        push!(G, sqrt((Approx[:,3][2])*sum((Approx[:,2]-Analytical[:,2]).^2)) )
    end
    
    return maximum(G)
    
end

function ComparisonPhaseRE(df,df2,T1,T2,T3,alpha)
    G=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)
        Astacked_df = [stack(ADFX[:,1:(end-1)])][1]
        Time=repeat(ADFX[:,end],Int(size(Astacked_df)[1]/size(ADFX)[1]))
        Astacked_df.Time = Time
        Analytical=Astacked_df
    
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)
        Bstacked_df = [stack(BDFX[:,1:(end-1)])][1]
        Time=repeat(BDFX[:,end],Int(size(Bstacked_df)[1]/size(BDFX)[1]))
        Bstacked_df.Time = Time
        Approx=Bstacked_df
        #push!(G,maximum(filter(!isnan,filter(!isinf,abs.(Approx[:,2]-Analytical[:,2])/abs.(Analytical[:,2])))))
        push!(G,maximum(abs.(Approx[:,2]-Analytical[:,2])./(1 .+ abs.(Analytical[:,2])) ))
    end 
    
    return maximum(G)
    
end

function ComparisonFrequencyRE(df,df2,T1,T2,T3,alpha)
    G=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)
        Astacked_df = [stack(ADFdX[:,1:(end-1)])][1]
        Time=repeat(ADFdX[:,end],Int(size(Astacked_df)[1]/size(ADFdX)[1]))
        Astacked_df.Time = Time
        Analytical=Astacked_df
    
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)
        Bstacked_df = [stack(BDFdX[:,1:(end-1)])][1]
        Time=repeat(BDFdX[:,end],Int(size(Bstacked_df)[1]/size(BDFdX)[1]))
        Bstacked_df.Time = Time
        Approx=Bstacked_df
        #push!(G,maximum(filter(!isnan,filter(!isinf,abs.(Approx[:,2]-Analytical[:,2])/abs.(Analytical[:,2])))))
        push!(G,maximum(abs.(Approx[:,2]-Analytical[:,2])./(1 .+ abs.(Analytical[:,2]))))
    end
    
    return maximum(G)
    
end

function ComparisonOH(df,df2,T1,T2,T3,alpha)
    G=[]
    H=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)

        ohi1=OverheatingIndicator(ADFX,df,df2)
    
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)

        ohi2=OverheatingIndicator(BDFX,df,df2)
        
        push!(G,ohi1)
        push!(H,ohi2)
    end 
    
    return maximum(G),maximum(H)
    
end


function ComparisonOHDF(df,df2,T1,T2,T3,alpha)
    G=[]
    H=[]
    STR=[]
    for i in 1:36
        Im=string(df.From[i])
        J=string(df.To[i])
        comma="#"
        str3="X"*Im*comma*J
        
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)

        ohi1=OverheatingIndicator(ADFX,df,df2)
    
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)

        ohi2=OverheatingIndicator(BDFX,df,df2)
        
        push!(G,ohi1)
        push!(H,ohi2)
        push!(STR,str3)
    end 

    A=fill(alpha,size(STR)[1])
    
    return G,H,STR,A
    
end 


function ComparisonOHREDEF(df,df2,T1,T2,T3,alpha)
    G=[]
    H=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)

        ohi1=OverheatingIndicatorREDEF(ADFX,df,df2)
    
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)

        ohi2=OverheatingIndicatorREDEF(BDFX,df,df2)
        
        push!(G,ohi1)
        push!(H,ohi2)
    end 
    
    return maximum(G),maximum(H)
    
end


function ComparisonOHDFREDEF(df,df2,T1,T2,T3,alpha)
    G=[]
    H=[]
    STR=[]
    for i in 1:36
        Im=string(df.From[i])
        J=string(df.To[i])
        comma="#"
        str3="X"*Im*comma*J
        
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)

        ohi1=OverheatingIndicatorREDEF(ADFX,df,df2)
    
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)

        ohi2=OverheatingIndicatorREDEF(BDFX,df,df2)
        
        push!(G,ohi1)
        push!(H,ohi2)
        push!(STR,str3)
    end 

    A=fill(alpha,size(STR)[1])
    
    return G,H,STR,A
    
end 




function ComparisonFlow(df,df2,T1,T2,T3,alpha)
    F=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)
        tdf=AnalysisGamma3(ADFX,df,df2)
        Ttdf=permutedims(tdf, 1)
        Ttdfstacked=[stack(Ttdf[:,3:end])][1]
        Time=repeat(Ttdf[:,2],Int(size(Ttdfstacked)[1]/size(Ttdf[:,3:end])[1]))
        Ttdfstacked.Time=Time
        rename!(Ttdfstacked, :value => "Anval")
        Analytical=Ttdfstacked
        tdf=AnalysisGamma3(BDFX,df,df2)
        Ttdf=permutedims(tdf, 1)
        Ttdfstacked=[stack(Ttdf[:,3:end])][1]
        Time=repeat(Ttdf[:,2],Int(size(Ttdfstacked)[1]/size(Ttdf[:,3:end])[1]))
        Ttdfstacked.Time=Time
        rename!(Ttdfstacked, :value => "Apval")
        Approx=Ttdfstacked
        Analytical.AbsFlow=abs.(Approx.Apval-Analytical.Anval)
        
        push!(F,maximum(Analytical.AbsFlow))
    end 

    return maximum(F)

end 

#=
function ComparisonFlowDF(df,df2,T1,T2,T3,alpha)
    L=[]
    AF=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)
        tdf=AnalysisGamma3(ADFX,df,df2)
        Ttdf=permutedims(tdf, 1)
        Ttdfstacked=[stack(Ttdf[:,3:end])][1]
        Time=repeat(Ttdf[:,2],Int(size(Ttdfstacked)[1]/size(Ttdf[:,3:end])[1]))
        Ttdfstacked.Time=Time
        rename!(Ttdfstacked, :value => "Anval")
        Analytical=Ttdfstacked
        tdf=AnalysisGamma3(BDFX,df,df2)
        Ttdf=permutedims(tdf, 1)
        Ttdfstacked=[stack(Ttdf[:,3:end])][1]
        Time=repeat(Ttdf[:,2],Int(size(Ttdfstacked)[1]/size(Ttdf[:,3:end])[1]))
        Ttdfstacked.Time=Time
        rename!(Ttdfstacked, :value => "Apval")
        Approx=Ttdfstacked
        Analytical.AbsFlow=abs.(Approx.Apval-Analytical.Anval)
        comp=combine(groupby(Analytical,:variable), :AbsFlow=>maximum=>:AbsFlow)
        line=comp[:,1]
        AbsFlow=comp[:,2]
        L=[L;line]
        AF=[AF;AbsFlow]
    end 

    A=fill(alpha,size(AF)[1])
    
    
    return L,AF,A
    
end 
=#

function ComparisonFlowDF(df,df2,T1,T2,T3,alpha)
    L=[]
    AF=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06)
        tdf=AnalysisGamma3(ADFX,df,df2)
        Ttdf=permutedims(tdf, 1)
        Ttdfstacked=[stack(Ttdf[:,3:end])][1]
        Time=repeat(Ttdf[:,2],Int(size(Ttdfstacked)[1]/size(Ttdf[:,3:end])[1]))
        Ttdfstacked.Time=Time
        rename!(Ttdfstacked, :value => "Anval")
        Analytical=Ttdfstacked
        tdf=AnalysisGamma3(BDFX,df,df2)
        Ttdf=permutedims(tdf, 1)
        Ttdfstacked=[stack(Ttdf[:,3:end])][1]
        Time=repeat(Ttdf[:,2],Int(size(Ttdfstacked)[1]/size(Ttdf[:,3:end])[1]))
        Ttdfstacked.Time=Time
        rename!(Ttdfstacked, :value => "Apval")
        Approx=Ttdfstacked
        Analytical.AbsFlow=abs.(Approx.Apval-Analytical.Anval)./(1 .+ abs.(Analytical.Anval) )
        comp=combine(groupby(Analytical,:variable), :AbsFlow=>maximum=>:AbsFlow)
        line=comp[:,1]
        AbsFlow=comp[:,2]
        L=[L;line]
        AF=[AF;AbsFlow]
    end 

    A=fill(alpha,size(AF)[1])
    
    
    return L,AF,A
    
end 

#=
T1=0
T2=0.5
T3=10


ALP=[0:0.01:1;]
L=[]
AF=[]
A=[]
for alpha in ALP
    #@time a=ComparisonPhase(df,df2,T1,T2,T3,alpha)
    #@time b=ComparisonFrequency(df,df2,T1,T2,T3,alpha)
    @time line,AbsFlow,A1=ComparisonFlowDF(df,df2,T1,T2,T3,alpha)
    L=[L;line]
    AF=[AF;AbsFlow]
    A=[A;A1]
end 

ABSComp=DataFrame(maxflow=AF,alpha=A,ind=L)
comp=combine(groupby(ABSComp,[:alpha,:ind]), :maxflow=>maximum=>:maxflow)
str="ABSComp10.csv"
CSV.write(str, comp)
=#


#=

https://medium.com/@mertbayraktarxd/what-is-smape-2cf605831feb

https://en.wikipedia.org/wiki/Symmetric_mean_absolute_percentage_error

https://math.stackexchange.com/questions/677852/how-to-calculate-relative-error-when-true-value-is-zero

ALP=[0:0.01:1;]
#alpha=0.90
#alpha=0.999999
T1=0
T2=0.5
T3=10

SMAPEComp=DataFrame(TrippedLine=[],maxSMAPE=[],alpha=[])


@time for alpha in ALP
    maxSMAPE=[]
    TrippedLine=[]
    for i in 1:36
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        ADFX=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.06)
        BDFX=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.06) 
        
        tdf=AnalysisGamma3(ADFX,df,df2)
        Ttdf=permutedims(tdf, 1)
        Ttdfstacked=[stack(Ttdf[:,3:end])][1]
        Time=repeat(Ttdf[:,2],Int(size(Ttdfstacked)[1]/size(Ttdf[:,3:end])[1]))
        Ttdfstacked.Time=Time
        rename!(Ttdfstacked, :value => "Anval")
        Analytical=Ttdfstacked
        tdf=AnalysisGamma3(BDFX,df,df2)
        Ttdf=permutedims(tdf, 1)
        Ttdfstacked=[stack(Ttdf[:,3:end])][1]
        Time=repeat(Ttdf[:,2],Int(size(Ttdfstacked)[1]/size(Ttdf[:,3:end])[1]))
        Ttdfstacked.Time=Time
        rename!(Ttdfstacked, :value => "Apval")
        Approx=Ttdfstacked
            
        Analytical.Apval=Approx.Apval      
        Analytical.AbsFlow=abs.(Approx.Apval-Analytical.Anval)
        Analytical.AbsFlowN=Analytical.AbsFlow./(abs.(Analytical.Anval)+abs.(Analytical.Apval))
        An=Analytical
        
        comp=combine(groupby(An,[:variable]), :AbsFlowN=>sum=>:AbsFlowN)
        comp.AbsFlowN=(comp.AbsFlowN)/(size(Ttdf[:,2])[1])
    
        push!(maxSMAPE,maximum(comp.AbsFlowN))
        TrippedLine=comp.variable
    
    end
    alp=fill(alpha,size(comp.variable)[1])
    sn=DataFrame(TrippedLine=comp.variable,maxSMAPE=maxSMAPE,alpha=alp)
    SMAPEComp=[SMAPEComp;sn]
end

=#






function EstimateGENOHI(Ntimes,df,df2,T1,T2,T3,noi=0.01,alpha=2/3,bs=0)
    i=rand(1:size(df)[1],1)[1]
    Df=copy(df)
    Df[i,4]=alpha*Df[i,4]
    Df2=copy(df2)
    DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
    #DFXS=AnalyticalSolution2(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
    oi=OverheatingIndicator(DFXS,df,df2)
    str="OverheatingIndicator"
    An=DataFrame(str=>oi)

    g=T1
    h=T2
    for j in 1:Ntimes
        i=rand(1:size(df)[1],1)[1]
        Df=copy(df)
        #Df[i,4]=0
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        # Duration of Failure is exponentailly distributed with lambda=T2-T1
        T2=g+rand(Exponential(h-T1),1)[1]
        DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
        #DFXS=AnalyticalSolution2(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
        oi=OverheatingIndicator(DFXS,df,df2)
        A=DataFrame(str=>oi)
        An=vcat(An,A)
    end
    return An 
end

function EstimateGENOHIindv3(df,df2,T1,T2,T3,noi=0.01,alpha=2/3,bs=0,delta=0.01)
    i=rand(1:size(df)[1],1)[1]
    Df=copy(df)
    Df[i,4]=alpha*Df[i,4]
    Df2=copy(df2)
    DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
    oi=OverheatingIndicator(DFXS,df,df2)
    str="OverheatingIndicator"
    df3=DataFrame(str=>oi)
    RT="T2"
    df4=DataFrame(RT=>T2)
    An=hcat(df3,df4)
    I=string(df.From[i])
    J=string(df.To[i])
    comma=","
    RL="X"*I*comma*J
    RemovedLine="RemovedLine"
    df5=DataFrame(RemovedLine=>RL)
    An=hcat(An,df5)
    for i=1:size(df)[1], T2=0:delta:1.5 #T2=0:0.001:1.5
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
        oi=OverheatingIndicator(DFXS,df,df2)
        df3=DataFrame(str=>oi)
        df4=DataFrame(RT=>T2)
        A=hcat(df3,df4)
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        RL="X"*I*comma*J
        df5=DataFrame(RemovedLine=>RL)
        A=hcat(A,df5)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end


function EstimateGENOHIindv(Ntimes,df,df2,T1,T2,T3,noi=0.01,alpha=2/3,bs=0)
    i=rand(1:size(df)[1],1)[1]
    Df=copy(df)
    Df[i,4]=alpha*Df[i,4]
    Df2=copy(df2)
    #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.0999)
    DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
    oi=OverheatingIndicatorIndv(DFXS,df,df2)
    
    #oi.T2=fill(T2,nrow(oi))
    #I=string(df.From[i])
    #J=string(df.To[i])
    #comma=","
    #RL="X"*I*comma*J
    #RemovedLine="RemovedLine"
    #oi.RemovedLine=fill(RL,nrow(oi))
    
    An=oi
    g=T1
    h=T2
    for j in 1:Ntimes
        i=rand(1:size(df)[1],1)[1]
        Df=copy(df)
        #Df[i,4]=0
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        # Duration of Failure is exponentailly distributed with lambda=T2-T1
        T2=g+rand(Exponential(h-T1),1)[1]
        #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.0999)
        DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
        oi=OverheatingIndicatorIndv(DFXS,df,df2)
        
        #oi.T2=fill(T2,nrow(oi))
        #I=string(df.From[i])
        #J=string(df.To[i])
        #comma=","
        #RL="X"*I*comma*J
        #RemovedLine="RemovedLine"
        #oi.RemovedLine=fill(RL,nrow(oi))
        A=oi
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end



function EstimateGENOHIindv2(df,df2,T1,T2,T3,noi=0.01,alpha=2/3,bs=0)
    i=rand(1:size(df)[1],1)[1]
    Df=copy(df)
    Df[i,4]=alpha*Df[i,4]
    Df2=copy(df2)
    DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
    oi=OverheatingIndicatorIndv(DFXS,df,df2)
    oi.T2=fill(T2,nrow(oi))
    I=string(df.From[i])
    J=string(df.To[i])
    comma=","
    RL="X"*I*comma*J
    RemovedLine="RemovedLine"
    oi.RemovedLine=fill(RL,nrow(oi))
    An=oi
    #for i=1:size(df)[1], T2=0:0.01:1.5
    for i=1:size(df)[1], T2=0.1:0.01:1.5
        Df=copy(df)
        Df[i,4]=alpha*Df[i,4]
        DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
        oi=OverheatingIndicatorIndv(DFXS,df,df2)
        oi.T2=fill(T2,nrow(oi))
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        RL="X"*I*comma*J
        RemovedLine="RemovedLine"
        oi.RemovedLine=fill(RL,nrow(oi))
        A=oi
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end










function EstimateGEN(Ntimes,df,df2,T1,T2,T3,noi=0.01,alpha=2/3,bs=0)
    i=rand(1:size(df)[1],1)[1]
    Df=copy(df)
    #alpha=0
    #alpha=2/3
    #Df[i,4]=0
    Df[i,4]=alpha*Df[i,4]
    Df2=copy(df2)
    #DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
    
    #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.09)
    #DFXS=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.09)
    oi=OverheatingIndicator(DFXS,df,df2)
    #oi=OverheatingIndicatorREDEF(DFXS,df,df2)
    str="OverheatingIndicator"
    df3=DataFrame(str=>oi)
    RT="T2"
    df4=DataFrame(RT=>T2)
    An=Analysis0(DFXS,df,df2)
    An=hcat(An,df3)
    An=hcat(An,df4)
    I=string(df.From[i])
    J=string(df.To[i])
    comma=","
    RL="X"*I*comma*J
    RemovedLine="RemovedLine"
    df5=DataFrame(RemovedLine=>RL)
    An=hcat(An,df5)
    g=T1
    h=T2
    for j in 1:Ntimes
        i=rand(1:size(df)[1],1)[1]
        Df=copy(df)
        #Df[i,4]=0
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        # Duration of Failure is exponentailly distributed with lambda=T2-T1
        T2=g+rand(Exponential(h-T1),1)[1]
        #DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.1,bs)
        
        #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.09)
        #DFXS=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.09)
        A=Analysis0(DFXS,df,df2)
        oi=OverheatingIndicator(DFXS,df,df2)
        #oi=OverheatingIndicatorREDEF(DFXS,df,df2)
        df3=DataFrame(str=>oi)
        df4=DataFrame(RT=>T2)
        A=hcat(A,df3)
        A=hcat(A,df4)
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        RL="X"*I*comma*J
        df5=DataFrame(RemovedLine=>RL)
        A=hcat(A,df5)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end




















function EstimateGENApprox(Ntimes,df,df2,T1,T2,T3,noi=0.01,alpha=2/3)
    i=rand(1:size(df)[1],1)[1]
    Df=copy(df)
    #alpha=0
    #alpha=2/3
    #Df[i,4]=0
    Df[i,4]=alpha*Df[i,4]
    Df2=copy(df2)
    #DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.07)
    
    #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.09)
    DFXS=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.099)
    oi=OverheatingIndicator(DFXS,df,df2)
    #oi=OverheatingIndicatorREDEF(DFXS,df,df2)
    str="OverheatingIndicator"
    df3=DataFrame(str=>oi)
    RT="T2"
    df4=DataFrame(RT=>T2)
    An=Analysis0(DFXS,df,df2)
    An=hcat(An,df3)
    An=hcat(An,df4)
    I=string(df.From[i])
    J=string(df.To[i])
    comma=","
    RL="X"*I*comma*J
    RemovedLine="RemovedLine"
    df5=DataFrame(RemovedLine=>RL)
    An=hcat(An,df5)
    g=T1
    h=T2
    for j in 1:Ntimes
        i=rand(1:size(df)[1],1)[1]
        Df=copy(df)
        #Df[i,4]=0
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        # Duration of Failure is exponentailly distributed with lambda=T2-T1
        T2=g+rand(Exponential(h-T1),1)[1]
        #DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.07)
        
        #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.09)
        DFXS=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.099)
        A=Analysis0(DFXS,df,df2)
        oi=OverheatingIndicator(DFXS,df,df2)
        #oi=OverheatingIndicatorREDEF(DFXS,df,df2)
        df3=DataFrame(str=>oi)
        df4=DataFrame(RT=>T2)
        A=hcat(A,df3)
        A=hcat(A,df4)
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        RL="X"*I*comma*J
        df5=DataFrame(RemovedLine=>RL)
        A=hcat(A,df5)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end



function EstimateGENnum(Ntimes,df,df2,T1,T2,T3,noi=0.01,alpha=2/3)
    i=rand(1:size(df)[1],1)[1]
    Df=copy(df)
    #alpha=0
    #alpha=2/3
    #Df[i,4]=0
    Df[i,4]=alpha*Df[i,4]
    Df2=copy(df2)
    #DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.02)
    
    #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.09)
    #DFXS=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.09)
    oi=OverheatingIndicator(DFXS,df,df2)
    #oi=OverheatingIndicatorREDEF(DFXS,df,df2)
    str="OverheatingIndicator"
    df3=DataFrame(str=>oi)
    RT="T2"
    df4=DataFrame(RT=>T2)
    An=Analysis0(DFXS,df,df2)
    An=hcat(An,df3)
    An=hcat(An,df4)
    I=string(df.From[i])
    J=string(df.To[i])
    comma=","
    RL="X"*I*comma*J
    RemovedLine="RemovedLine"
    df5=DataFrame(RemovedLine=>RL)
    An=hcat(An,df5)
    g=T1
    h=T2
    for j in 1:Ntimes
        i=rand(1:size(df)[1],1)[1]
        Df=copy(df)
        #Df[i,4]=0
        Df[i,4]=alpha*Df[i,4]
        Df2=copy(df2)
        # Duration of Failure is exponentailly distributed with lambda=T2-T1
        T2=g+rand(Exponential(h-T1),1)[1]
        #DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.02)
        
        #DFXS=AnalyticalSolution(df,df2,Df,Df2,T1,T2,T3,0.09)
        #DFXS=AnalyticalSolutionApprox(df,df2,Df,Df2,T1,T2,T3,0.09)
        A=Analysis0(DFXS,df,df2)
        oi=OverheatingIndicator(DFXS,df,df2)
        #oi=OverheatingIndicatorREDEF(DFXS,df,df2)
        df3=DataFrame(str=>oi)
        df4=DataFrame(RT=>T2)
        A=hcat(A,df3)
        A=hcat(A,df4)
        I=string(df.From[i])
        J=string(df.To[i])
        comma=","
        RL="X"*I*comma*J
        df5=DataFrame(RemovedLine=>RL)
        A=hcat(A,df5)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end





function Estimate0(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    oi=OverheatingIndicator(DFXS,df,df2)
    str="OverheatingIndicator"
    df3=DataFrame(str=>oi)
    RT="T2"
    df4=DataFrame(RT=>T2)
    An=Analysis0(DFXS,df,df2)
    An=hcat(An,df3)
    An=hcat(An,df4)
    g=T1
    h=T2
    for i in 1:Ntimes
        # Duration of Failure is exponentailly distributed with lambda=T2-T1
        T2=g+rand(Exponential(h-T1),1)[1]
        DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        A=Analysis0(DFXS,df,df2)
        oi=OverheatingIndicator(DFXS,df,df2)
        df3=DataFrame(str=>oi)
        df4=DataFrame(RT=>T2)
        A=hcat(A,df4)
        A=hcat(A,df3)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end




function Estimate(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    An=Analysis(DFXS,df,df2)
    h=T2
    for i in 1:Ntimes
        T2=h+rand(Exponential(h-T1),1)[1]
        DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        A=Analysis(DFXS,df,df2)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end

function EstimateBeta(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    An=AnalysisBeta(DFXS,df,df2)
    h=T2
    for i in 1:Ntimes
        T2=h+rand(Exponential(h-T1),1)[1]
        DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        A=AnalysisBeta(DFXS,df,df2)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end

function EstimateDelta(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    An=AnalysisDelta(DFXS,df,df2)
    h=T2
    for i in 1:Ntimes
        T2=h+rand(Exponential(h-T1),1)[1]
        DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        A=AnalysisDelta(DFXS,df,df2)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end

function EstimateEpsilon(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    An=AnalysisEpsilon(DFXS,df,df2)
    h=T2
    for i in 1:Ntimes
        T2=h+rand(Exponential(h-T1),1)[1]
        DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        A=AnalysisEpsilon(DFXS,df,df2)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end

function EstimateZeta(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
    An=AnalysisZeta(DFXS,df,df2)
    h=T2
    for i in 1:Ntimes
        T2=h+rand(Exponential(h-T1),1)[1]
        DFXS=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        A=AnalysisZeta(DFXS,df,df2)
        An=vcat(An,A)
    end
    #CSV.write("An.csv",An)
    return An 
end

function Analysis2(Kappa,df,df2)
    DF=filter(:Simulation => isequal(Kappa.Simulation[1]), Kappa)
    DF=select!(DF, Not(:Simulation))
    An=Analysis(DF,df,df2)[1]
    for i in unique(Kappa.Simulation)
        D=filter(:Simulation => isequal(i), Kappa)
        D=select!(D, Not(:Simulation))
        A=Analysis(D,df,df2)
        An=vcat(An,A)
        
    end
    return An
end

function Estimate2(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    Simulation="Simulation"
    DFXS=SimulationDEFSDE(df,df2,Df,Df2,T1,T2,T3,noi)
    df3=DataFrame(Simulation=>fill("0", nrow(DFXS)))
    DFXS=hcat(DFXS,df3)
    for i in 1:Ntimes
        I=string(i)
        Simulation="Simulation"
        D=SimulationDEFSDE2(df,df2,Df,Df2,T1,T2,T3,noi)
        df3=DataFrame(Simulation=>fill(I, nrow(D)))
        D=hcat(D,df3)
        DFXS=vcat(DFXS,D)
    end
    CSV.write("DFXS.csv",DFXS)
    return DFXS 
end



# A Fault in the System (Lines) + Clearance ODE
function SimulationDEFODE(df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    n=df.Nnodes[1] #Number of of buses    
    P1=df2.PowerInjections #Vector of power injections 1
    A3,A4=XiMass(df,df2) #Xi and Mass Matrix 1
    function f(du,u,p,t)
        du.=A4*u+[P1;zeros(n)]
    end
    L1=Laplatian(df) # Laplatian matrix 1
    Theta0=pinv(L1)*P1 #Steady-state solution
    
    u0 = [zeros(n);Theta0];
    #u0 = [Theta0;zeros(n)];
    
    tspan1 = (0.0,T1);
    ###
    #ODE 1
    fun = ODEFunction(f; mass_matrix=A3);
    prob1 = ODEProblem(fun,u0,tspan1);
    sol1 = solve(prob1,Rosenbrock23());
    ###
    RES = sol1.u[1];
    for i in 1:length(sol1.t)
       RES=hcat(RES,sol1.u[i]) 
    end
    P2=Df2.PowerInjections #Vector of power injections 2
    A3,A4=XiMass(Df,Df2) #Xi and Mass Matrix 2
    function f2(du,u,p,t)
        du.=A4*u+[P2;zeros(n)]
    end
    # Fault in the system
    u2 =sol1.u[end]
    tspan2 = (T1,T2);
    ###
    #ODE 2
    fun2 = ODEFunction(f2; mass_matrix=A3);
    prob2 = ODEProblem(fun2,u2,tspan2);
    sol2 = solve(prob2,Rosenbrock23());
    ###
    for i in 1:length(sol2.t)
       RES=hcat(RES,sol2.u[i]) 
    end
    #Clearance
    u3 =sol2.u[end]
    tspan3 = (T2,T3);
    A3,A4=XiMass(df,df2);
    function f3(du,u,p,t)
        du.=A4*u+[P1;zeros(n)]
    end
    fun3 = ODEFunction(f3; mass_matrix=A3);
    ###
    #ODE 3
    prob3 = ODEProblem(fun3,u3,tspan3);
    sol3 = solve(prob3,Rosenbrock23());
    ###
    for i in 1:length(sol3.t)
       RES=hcat(RES,sol3.u[i]) 
    end
    #Store Info as DataFrame
    MAT=RES[:,2:end]
    J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    #DFdX=DataFrame(J[:,1:n],:auto)
    #=
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:n#size(DFX)[2]
        for j in 1:n#size(DFX)[2]
            if j>i
                I=string(i)
                J=string(j)
                comma=","
                str="X"*I*comma*J
                df3=DataFrame(str=>DFX[:,i]-DFX[:,j])
                DFX=hcat(DFX,df3)
                                                
            end
        end
    end
    =#
    DFX.Time = range(0, T3, length=size(DFX)[1])
    #CSV.write("DFX1.csv", DFX)
    #CSV.write("DFdX.csv", DFdX)
    return DFX
end





function SimulationDEFODE2(df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    
    n=df.Nnodes[1] #Number of of buses    
    P1=df2.PowerInjections #Vector of power injections 1
    A3,A4=XiMass(df,df2) #Xi and Mass Matrix 1 (Initial Confiduration)
    P2=Df2.PowerInjections #Vector of power injections 2
    A5,A6=XiMass(Df,Df2) #Xi and Mass Matrix 2 (Fault)
    function f(du,u,p,t)
        if (t>=T1)&(t<=T2)
            du.=A6*u+[P2;zeros(n)]
        else
            du.=A4*u+[P1;zeros(n)]
        end
    end
    L1=Laplatian(df) # Laplatian matrix 1
    Theta0=pinv(L1)*P1 #Steady-state solution
    
    u0 = [zeros(n);Theta0];
    #u0 = [Theta0;zeros(n)];
    
    tspan1 = (0.0,T3);
    ###
    #ODE 
    fun = ODEFunction(f; mass_matrix=A3);
    #SDEProblem(SDEFunction(f,noise; mass_matrix=A3),u0,tspan1,noise)
    prob1 = ODEProblem(fun,u0,tspan1);
    sol1 = solve(prob1,Rodas5P());
    ###
    RES = sol1.u[1];
    for i in 1:length(sol1.t)
       RES=hcat(RES,sol1.u[i]) 
    end
    #Store Info as DataFrame
    MAT=RES[:,2:end]
    J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    #DFdX=DataFrame(J[:,1:n],:auto)
    #=
    #Compute X_{ij}=X_i-X_j for all i,j=1,...,n
    for i in 1:n#size(DFX)[2]
        for j in 1:n#size(DFX)[2]
            if j>i
                I=string(i)
                J=string(j)
                comma=","
                str="X"*I*comma*J
                df3=DataFrame(str=>DFX[:,i]-DFX[:,j])
                DFX=hcat(DFX,df3)
                                                
            end
        end
    end
    =#
    DFX.Time = range(0, T3, length=size(DFX)[1])
    #CSV.write("DFX1.csv", DFX)
    #CSV.write("DFdX.csv", DFdX)
    return DFX
end

# A Fault in the System (Lines) + Clearance SDE
function SimulationDEFSDE(df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    #T0=0 init
    #T1 fault time
    #T2 clearance
    #T3 finish
    n=df.Nnodes[1] #Number of of buses    
    P1=df2.PowerInjections #Vector of power injections 1
    A3,A4=XiMass(df,df2) #Xi and Mass Matrix 1
    function f(du,u,p,t)
        du.=A4*u+[P1;zeros(n)]
    end
    #Noise for SDE
    function noise(du, u, p, t)
        for i in 1:n
            du[i] = noi
        end 
    end
    L1=Laplatian(df) # Laplatian matrix 1
    Theta0=pinv(L1)*P1 #Steady-state solution
    u0 = [zeros(n);Theta0];
    tspan1 = (0.0,T1);
    #SDE 1
    prob1 = SDEProblem(SDEFunction(f,noise; mass_matrix=A3),u0,tspan1,noise);
    sol1 = solve(prob1,ImplicitEM());
    RES = sol1.u[1];
    for i in 1:length(sol1.t)
       RES=hcat(RES,sol1.u[i]) 
    end
    P2=Df2.PowerInjections #Vector of power injections 2
    A3,A4=XiMass(Df,Df2) #Xi and Mass Matrix 2
    function f2(du,u,p,t)
        du.=A4*u+[P2;zeros(n)]
    end
    # Fault in the system
    u2 =sol1.u[end]
    tspan2 = (T1,T2);
    #SDE 2
    prob2 = SDEProblem(SDEFunction(f2,noise; mass_matrix=A3),u2,tspan2,noise);
    sol2 = solve(prob2,ImplicitEM());
    ###
    for i in 1:length(sol2.t)
       RES=hcat(RES,sol2.u[i]) 
    end
    #Clearance
    u3 =sol2.u[end]
    tspan3 = (T2,T3);
    A3,A4=XiMass(df,df2) #Xi and Mass Matrix 1
    function f3(du,u,p,t)
        du.=A4*u+[P1;zeros(n)]
    end
    #SDE 3
            #SDEProblem(SDEFunction(f,noise; mass_matrix=A3),u0,tspan1,noise)
    prob3 = SDEProblem(SDEFunction(f3,noise; mass_matrix=A3),u3,tspan3,noise);
    sol3 = solve(prob3,ImplicitEM());
    ###
    for i in 1:length(sol3.t)
       RES=hcat(RES,sol3.u[i]) 
    end
    #Store Info as DataFrame
    MAT=RES[:,2:end]
    J=transpose(MAT)
    DFX=DataFrame(J[:,n+1:2*n],:auto)
    #DFdX=DataFrame(J[:,1:n],:auto)
    #CSV.write("DFX2.csv", DFX)
    #CSV.write("DFdX.csv", DFdX)
    DFX.Time = range(0, T3, length=size(DFX)[1])
    return DFX
end




















#=
function PLT1(DF,v1,v2,v3,X,Y,tspan)
    @manipulate for 
    b12 in slider(0.01:0.01:10.0, value=10, show_value=true,label="b12"),
    b13 in slider(0.01:0.01:10.0, value=10, show_value=true,label="b13"),
    b23 in slider(0.01:0.01:10.0, value=10, show_value=true,label="b23"),
    f12 in slider(0.01:0.01:10.0, value=v1, show_value=true,label="f12 max"),
    f13 in slider(0.01:0.01:10.0, value=v2, show_value=true,label="f13 max"),
    f23 in slider(0.01:0.01:10.0, value=v3, show_value=true,label="f23 max")
    th12max = f12/b12 
    th13max = f13/b13
    th23max = f23/b23 
    Mean1 = mean(Bool.(abs.(DF[:, "X1,2"]).<=th12max).&&
                 Bool.(abs.(DF[:, "X1,3"]).<=th13max).&&
                 Bool.(abs.(DF[:, "X2,3"]).<=th23max))
    Mean2 = tspan[2]*Mean1
    Mean3 = tspan[2] - Mean2
    ARG = Bool.(abs.(DF[:, "X1,2"]).>th12max).||
        Bool.(abs.(DF[:, "X1,3"]).>th13max).||
        Bool.(abs.(DF[:, "X2,3"]).>th23max)
    if sum(ARG)>0
            Time1 = minimum(DF.Time[ARG])
    else
            Time1="Undefined"
    end
    B = [b12+b13 -b12 -b13; -b12 b13+b23 -b23; -b13 -b23 b13+b23]
    plot([th12max;th12max],[-pi;pi])
    plot!([-th12max;-th12max],[-pi;pi])
    plot!([-pi;pi],[th13max;th13max])
    plot!([-pi;pi], [-th13max;-th13max])
    plot!([-pi;pi],[th23max-pi;th23max+pi])
    plot!([-pi;pi], [-th23max-pi;-th23max+pi])
    plot!(legend=false, xlabel="θ12", ylabel="θ13",
            xlim=X,ylim=Y)
    scatter!(DF[:, "X1,2"],DF[:, "X1,3"])
    annotate!(0, 0.5, text("Probability of being inside the polytope: $Mean1", :green, :center, 8))
    annotate!(0, 0.4, text("Mean time inside the polytope: $Mean2", :blue, :center, 8))    
    annotate!(0, 0.3, text("Mean time outside the polytope: $Mean3", :red, :center, 8))
    annotate!(0, 0.2, text("Fist time outside the polytope: $Time1", :purple, :center, 8))    
    end
end

function PLT2(DF)
    scatter!(DF[:, "X1,2"],DF[:, "X1,3"])
end

function PLT3(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    DFXS=SimulationDEFSDE(df,df2,Df,Df2,T1,T2,T3,noi)
    for i in 1:Ntimes
        D=SimulationDEFSDE(df,df2,Df,Df2,T1,T2,T3,noi)
        DFXS=vcat(DFXS,D)
    end
    plot(DFXS[:, "X1,2"],DFXS[:, "X1,3"])
end 

function PLT4(Ntimes,df,df2,Df,Df2,T1,T2,T3,noi=0.01)
    DFXS=SimulationDEFSDE(df,df2,Df,Df2,T1,T2,T3,noi)
    for i in 1:Ntimes
        D=SimulationDEFSDE(df,df2,Df,Df2,T1,T2,T3,noi)
        DFXS=vcat(DFXS,D)
    end
    scatter(DFXS[:, "X1,2"],DFXS[:, "X1,3"])
end 
=#