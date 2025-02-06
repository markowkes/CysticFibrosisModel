using Interpolations
using CSV, DataFrames

include("solver.jl")

p = (
    Name = "CaseA",
    μ   = 0.1,         # Growth rate in mucus layer [1/h]
    b   = 0.2,         # Death/decay rate in mucus layer [1/h]
    v   = 6,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 0,           # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> 10^4, # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 20,       # Final time
)
solA = solver(p)

p = (
    Name = "CaseB",
    μ   = 0.1,         # Growth rate in mucus layer [1/h]
    b   = 0.2,         # Death/decay rate in mucus layer [1/h]
    v   = 0,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 0,           # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> 10^4, # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 100,      # Final time
)
solB = solver(p)

p = (
    Name = "CaseC",
    μ   = 0.3,         # Growth rate in mucus layer [1/h]
    b   = 0,           # Death/decay rate in mucus layer [1/h]
    v   = 6,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 0,           # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> 10^4, # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 100,      # Final time
)
solC = solver(p)

p = (
    Name = "CaseD",
    μ   = 0.3,         # Growth rate in mucus layer [1/h]
    b   = 0,           # Death/decay rate in mucus layer [1/h]
    v   = 0,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 0,           # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> 10^4, # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 100,      # Final time
)
solD = solver(p)

p = (
    Name = "CaseE",
    μ   = 0.3,         # Growth rate in mucus layer [1/h]
    b   = 0,           # Death/decay rate in mucus layer [1/h]
    v   = 3,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 10^4,        # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> 10^4, # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 1000,      # Final time
)
solE = solver(p)

# Get steady state solution from E 
z,dz = createGrid(p)
SSE = [p.Xo; solE(p.tFinal)]
interp=linear_interpolation(z,SSE)

p = (
    Name = "CaseF",
    μ   = 0.3,         # Growth rate in mucus layer [1/h]
    b   = 0,           # Death/decay rate in mucus layer [1/h]
    v   = 3,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 10^4,        # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> 0,    # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 100,      # Final time
)
solF = solver(p)

p = (
    Name = "CaseG",
    μ   = 0.1,         # Growth rate in mucus layer [1/h]
    b   = 0.2,           # Death/decay rate in mucus layer [1/h]
    v   = 6,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 10^4,        # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> interp(z), # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 100,      # Final time
)
solG = solver(p)

p = (
    Name = "CaseH",
    μ   = 0.1,         # Growth rate in mucus layer [1/h]
    b   = 0.2,           # Death/decay rate in mucus layer [1/h]
    v   = 6,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 10^4,        # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> 0,    # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 100,      # Final time
)
solH = solver(p)

p = (
    Name = "CaseI",
    μ   = 0.3,         # Growth rate in mucus layer [1/h]
    b   = 0,           # Death/decay rate in mucus layer [1/h]
    v   = 3,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 10^3,        # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> interp(z) , # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 100,      # Final time
)
solI = solver(p)

p = (
    Name = "CaseJ",
    μ   = 0.1,         # Growth rate in mucus layer [1/h]
    b   = 0.2,           # Death/decay rate in mucus layer [1/h]
    v   = 6,           # Velocity of mucocilliary escalator [cm/h]
    Xo  = 10^3,        # Concentration at bottom of lung (z=0) [cfu/cm³]
    Xin = (z) -> interp(z) , # Initial concentration [cfu/cm³]
    L   = 48,          # Length of the escalator [cm]
    Nz  = 1000,          # Number of grid points
    tFinal = 100,      # Final time
)
solJ = solver(p)
