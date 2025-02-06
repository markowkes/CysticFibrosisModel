using DifferentialEquations
using Plots
using LaTeXStrings
using Printf

function solver(p)
    println("Starting solver: ",p.Name)

    z,dz = createGrid(p)

    # Prepare IC 
    X0 = zeros(p.Nz)
    for i in 1:p.Nz
        X0[i] = p.Xin(z[i+1])
    end

    # Prepare ODE solver 
    prob = ODEProblem(RHS!,X0,(0.0,p.tFinal),p)

    sol = solve(prob)

    processSol(sol,p)

    return sol
end

function RHS!(dsol,sol,p,t)
    z,dz = createGrid(p)

    # Compute RHS 
    for i in 1:p.Nz 
        Xm = i==1 ? p.Xo : sol[i-1]
        dXdz = (sol[i] - Xm)/dz
        dsol[i] = (p.μ - p.b)*sol[i] - p.v*dXdz
    end
end

function createGrid(p)
    z = LinRange(0,p.L,p.Nz+1)
    dz = z[2]-z[1]
    return z,dz
end

function processSol(sol,p)
    plot_Xtop_vs_Time(sol,p)
    plot_X_vs_z_movie(sol,p)
    sol2CSV(sol,p)
end

function plot_Xtop_vs_Time(sol,p)
    # Number of times to plot
    N = 100
    # Create times
    times = LinRange(sol.t[1],sol.t[end],N)
    # Preallocate X
    Xtop = zeros(N)
    # Get X at top of lung at each time
    for i in 1:N
        Xtop[i] = sol(times[i])[end]
    end
    plt = plot(times,Xtop, legend=false)
    xaxis!(L"\textrm{Time}~[\mathrm{h}]")
    yaxis!(L"\textrm{Top~Concentration}~X_{L}(t)~[\mathrm{cfu}~\mathrm{cm}^{-3}]")
    ylims!(0,1.1*maximum(maximum(sol)))
    display(plt)
    savefig(p.Name*".png")
end

function plot_X_vs_z(sol,t,p)
    z,dz = createGrid(p)
    X = [p.Xo;sol(t)]
    plt = plot(z,X, legend=false)
    xaxis!(L"\textrm{Position}~[\mathrm{cm}]")
    yaxis!(L"\textrm{Concentration}~X(z)~[\mathrm{cfu}~\mathrm{cm}^{-3}]")
    plot_title = @sprintf("t = %.2f",t)
    ylims!(0,maximum(maximum(sol)))
    title!(plot_title)
    display(plt)
end

function plot_X_vs_z_movie(sol,p)
    # Number of times to plot
    N = 100
    # Create times
    times = LinRange(sol.t[1],sol.t[end],N)
    anim = @animate for i in eachindex(times)
        plot_X_vs_z(sol,times[i],p)
    end
    gif(anim, p.Name*".gif", fps=15)
end

function sol2CSV(sol,p)
    # Number of times to write
    N = 1000
    # Create times
    times = LinRange(sol.t[1],sol.t[end],N)
    # Preallocate X
    Xtop = zeros(N)
    # Get X at top of lung at each time
    for i in 1:N
        Xtop[i] = sol(times[i])[end]
    end
    
    data = DataFrame(time = times,
                     Xtop = Xtop)
    CSV.write(p.Name*".csv",data)
end
