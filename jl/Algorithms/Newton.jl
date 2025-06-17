using NonlinearSolve
using DynamicalSystems
using ForwardDiff
using DifferentialEquations
using LinearAlgebra
using StaticArrays


function newton_method(ds::CoupledODEs, u0::U where {U<:AbstractArray{<:Real}}, T::Float64, Δt::Float64, t_min::Float64)
    # Dimensionality of the system
    D = dimension(ds)

    # Let's solve the problem for systems of dimension 3
    @assert(D == 3)

    # v - the value of the phase velocity vector at point u0
    v = ds.integ.f.f(u0, ds.p0, 0.0)

    # k - index of the component that determines the position of the secant
    k = argmax(abs.(v))

    # i and j - indices of the remaining components that will change
    i, j = filter((x) -> (return x != k), [i for i in 1:D])[1:2]

    # Estimated time of the period
    T_local = ForwardDiff.value(x[end])

    # Ranges in which we will integrate the system
    tspan = (0.0, T_local)

    # Detefinition of the residual
    f = (err, x, p) -> begin
      # Creating a vector of initial conditions
      u = SVector{D}(
        if k == 1
          (u0[k], x[1], x[2])
        elseif k == 2
          (x[1], u0[k], x[2])
        elseif k == 3
          (x[1], x[2], u0[k])
        end
      )

      # Integration of the system for the initial condition u for time T_local
      sol = solve(remake(ds.integ.sol.prob; u0 = u, tspan = tspan);
                  DynamicalSystemsBase.DEFAULT_DIFFEQ..., ds.diffeq...)

      # The difference between the final and initial state
      err = sol.u[end] .- sol.u[1]
    end

    # Inital state
    x0 = [u0[[i, j]]; T]

    prob = NonlinearProblem(f, x0, nothing)

    sol = solve(prob, NewtonRaphson();
                show_trace = Val(true), trace_level = TraceAll())

    u0_new = zeros(3)
    u0_new[i] = sol.u[1]
    u0_new[j] = sol.u[2]
    u0_new[k] = u0[k]

    T_new  = sol.u[end]

    T_new < t_min && (println("Newton: Period $T_new is less than t_min"); return nothing)

    return (state = trajectory(ds, T_new - Δt, u0_new; Δt = Δt)[1], period = T_new)
end

function main()
  system = Systems.lorenz()

  ig = [-1.9598067770796535, 2.0076976037190195, 26.58579174986873]
  T  = 39.44
  Δt = 0.01
  t_min = 1.0

  res = trajectory(system, T - Δt, ig; Δt = Δt)[1]
  Δ   = LinearAlgebra.norm(res[end] .- res[1])

  res_newton = newton_method(system, ig, T, Δt, t_min)
  !isnothing(res_newton) && (Δ_newton = LinearAlgebra.norm(res_newton.state[end] .- res_newton.state[1]))

  res_lm = lm(system, ig, T, Δt, t_min)
  !isnothing(res_lm) && (Δ_lm = LinearAlgebra.norm(res_lm.state[end] .- res_lm.state[1]))

  println("Begin state = $ig")
  println("Δ = $Δ, T = $T")
  println("res[end] - res[1] = $(res[end] .- res[1])")
  !isnothing(res_newton) && (println("Newton begin state = $(res_newton.state[1])"))
  !isnothing(res_newton) && (println("Δ_newton = $Δ_newton, T = $(res_newton.period)"))
  !isnothing(res_lm) && (println("LM begin state = $(res_lm.state[1])"))
  !isnothing(res_lm) && (println("Δ_lm = $Δ_lm, T = $(res_lm.period)"))
end

main()
