#using Pkg; Pkg.add(url = "https://github.com/JuliaDynamics/PeriodicOrbits.jl")
using PeriodicOrbits
using NonlinearSolve
using SciMLBase


"""
This function is almost full copy of PeriodicOrbits.periodic_orbit()
Reason why it's here is because the original one doesn't support Δt changes
while creating PO trajectory.
"""
function periodic_orbit_(ds::CoupledODEs, alg::OptimizedShooting, ig::InitialGuess, Δt::Float64, t_min::Float64)
  D = dimension(ds)

  f = (err, v, p) -> begin
    if isinplace(ds) 
        u0 = @view v[1:D]
    else
        u0 =  SVector{D}(v[1:D])
    end
    T = v[end]

    bounds = zeros(eltype(v), alg.n * 2)
    for i in 0:alg.n-1
        bounds[i + 1] = i * alg.Δt
        bounds[i + alg.n + 1] = T + i * alg.Δt
    end
    tspan = (0.0, T + alg.n * alg.Δt)

    sol = solve(SciMLBase.remake(ds.integ.sol.prob; u0 = u0, tspan = tspan);
                DynamicalSystemsBase.DEFAULT_DIFFEQ..., ds.diffeq..., saveat = bounds)
    if (length(sol.u) == alg.n * 2)
        for i in 1:alg.n
            err[D * i - (D - 1):D * i] = (sol.u[i] - sol.u[i + alg.n])
        end
    else
        fill!(err, Inf)
    end
  end

  prob = NonlinearLeastSquaresProblem(
      NonlinearFunction(f, resid_prototype = zeros(alg.n*dimension(ds))), [ig.u0..., ig.T])

  sol = solve(prob, NonlinearSolve.LevenbergMarquardt(); alg.nonlinear_solve_kwargs...)

  u0 = sol.u[1:end - 1]
  T  = sol.u[end]

  T < t_min && (return nothing)

  return (state = trajectory(ds, T - Δt, u0; Δt = Δt)[1], period = T)
end


function lm(system::CoupledODEs, u0::U where {U<:AbstractArray{<:Real}}, T::Float64, Δt::Float64, t_min::Float64)
  ig = InitialGuess(u0, T)
  po = periodic_orbit_(system, OptimizedShooting(), ig, Δt, t_min)
  return po
end
