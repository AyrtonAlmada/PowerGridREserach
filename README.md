# Real-Time Stochastic Assessment of Dynamic N-1 Grid Contingencies

This repository contains a collection of Julia implementations for simulating and analyzing power grid dynamics under transmission line faults. The framework is designed to study the transient behavior of power systems, including the evolution of generator phase angles and frequency deviations following disturbances such as line outages.

The code implements both analytical solutions of linearized swing dynamics and stochastic simulations, allowing researchers to investigate how faults propagate through the network and how long transmission lines operate beyond their thermal limits.

The toolkit was developed as part of research on grid stability, cascading failures, and overheating risk in transmission networks.

## Overview

Power transmission networks are modeled as weighted graphs where nodes represent buses (generators or loads) and edges represent transmission lines. The dynamic behavior of the system is described using the classical swing equation, which captures the interaction between generator inertia, damping, and network coupling.

The repository provides tools for:

* constructing the network Laplacian from transmission data

* building mass and system matrices of the swing equation

* computing analytical time-domain solutions during fault events

* performing stochastic simulations of grid dynamics

* evaluating overheating indicators for transmission lines

* detecting constraint violations and fault propagation

* conducting Monte-Carlo experiments to estimate system risk

These tools enable the study of transient stability and thermal stress in power grids following line failures.

## Key Features

### Network Modeling

* Construction of the weighted Laplacian matrix representing the transmission grid topology.

* Support for generator parameters including inertia, damping, and power injections.

### Dynamic Simulation

Two types of simulations are available:

#### Analytical simulation

* Closed-form solution of the linearized swing dynamics.

* Efficient evaluation of trajectories during fault and post-fault phases.

#### Stochastic simulation

* Euler–Maruyama integration of stochastic swing dynamics.

* Random perturbations representing load fluctuations or renewable variability.

### Fault Modeling

The framework models line faults and clearance events using a three-stage timeline:

```
t < T1      : normal operation
T1 ≤ t ≤ T2 : faulted network topology
t > T2      : post-fault recovery
```

The transmission network is modified during the fault interval to represent the removal or degradation of a line.

### Overheating Indicators

The repository implements metrics that quantify thermal stress in transmission lines by measuring the time during which phase differences exceed safe limits.

These indicators provide a proxy for line overheating risk and are useful for studying cascading failures.

### Monte-Carlo Analysis

Random fault scenarios can be generated to estimate the distribution of overheating indicators and identify vulnerable lines in the network.

## Repository Structure

```
PowerGridDynamics/
│
├── 0src/
│   └── PowerGridsFunctions3.jl
│
├── 0examples/
│   ├── SimulationExampleISRAEL.ipynb
│   ├── IsraeliPowerGrid.png
│   └── israel.m
│
├── 0data/
│   ├── branch.csv
│   ├── branch2.csv
│   └── bus.csv
│
├── 1README/
│   ├── 0AyrtonAlmadaResume.pdf
│   ├── 1AyrtonAlmadaPaper.pdf
│   └── 2AyrtonAlmadaPoster.pdf
│
├── 1SingPhase/
│   └── [Examples single phase faults (partial removal) in the Israeli system]
│       └── [[For each case we show Evolution of Frequencies, Normalized Power Flows, Phases,
│       └─────Overview Analysis (mp4/gif), Overloaded Lines of the System and Overheating. ]]
│
├── 1ThreePhases/
│   └── [Examples three phase faults (complete removal) in the Israeli system]
│       └── [[For each case we show Evolution of Frequencies, Normalized Power Flows, Phases,
│       └─────Overview Analysis (mp4/gif), Overloaded Lines of the System and Overheating. ]]
│
└── README.md
```

## Installing

Ensure that Julia (version ≥ 1.8) is installed.

Clone the repository:

```
git clone https://github.com/yourusername/power-grid-dynamics.git
cd power-grid-dynamics
```

Install dependencies in Julia:

```
using Pkg
Pkg.add([
    "DataFrames",
    "LinearAlgebra",
    "Statistics",
    "SparseArrays",
    "Distributions",
    "Roots"
])
```

## Example Usage

```
branch=DataFrame(CSV.File("branch.csv"));
bus=DataFrame(CSV.File("bus.csv"));
# Re-label data
label2id=Dict{Int,Int}()
for (i,j) in enumerate(bus.index)
    push!(label2id,j => i)
end
indexf_bus=[label2id[i] for i in branch.f_bus]
indext_bus=[label2id[i] for i in branch.t_bus]
indexbus=[label2id[i] for i in bus.index];
bus.indexbus=indexbus
branch.indexf_bus=indexf_bus
branch.indext_bus=indext_bus;
# Drop orginal labels
select!(branch, Not([:f_bus, :t_bus, :index]));
df=DataFrame(Lines = (1:nrow(branch)), From = branch.indexf_bus, To = branch.indext_bus, 
    Susceptance = (1)./(branch.br_x), F = branch.rate_a, Nnodes = fill(nrow(bus), nrow(branch)))
df.ThetaMax=df.F./df.Susceptance
df=unique(df,[:From,:To])
df=sort!(df,[:From,:To])
df.Lines=(1:nrow(df))
df2=DataFrame(Damping = bus.d, Inertia = bus.m, PowerInjections = bus.p, Theta0 = bus.va);

T1=0.0 #Time of Failure
T2=5 #Time of Clearance
T3=50 #Finish time
Bf=copy(df)
j=1
Bf[j,4]=0*Bf[j,4]
Bf2=copy(df2)
@time DFX=AnalyticalSolution(df,df2,Bf,Bf2,T1,T2,T3,0.01)
@time b=OverheatingIndicator(DFX,df,df2)
println("Overheating indicator: ", b)
```

## Applications

This toolkit can be used for research in:

* power system transient stability

* grid resilience and cascading failures

* renewable integration and stochastic disturbances

* thermal stress analysis of transmission networks

* risk assessment under random faults

## Citation

If you use this code in your research, please cite the following publications:

```
@article{almada2025real,
  title={Real-Time Stochastic Assessment of Dynamic N-1 Grid Contingencies},
  author={Almada, Ayrton and Pagnier, Laurent and Goldshtein, Igal and Kazi, Saif R and others},
  journal={arXiv preprint arXiv:2510.18007},
  year={2025}
}
```

```
@article{almada2026real,
  title={Real-Time Dynamic N-1 Screening: Identifying High-Risk Lines and Transformers After Common Faults},
  author={Almada, Ayrton and Pagnier, Laurent and Goldshtein, Igal and Kazi, Saif R and others},
  journal={arXiv preprint arXiv:2602.12293},
  year={2026}
}
```

## License

This project is released under the MIT License.








