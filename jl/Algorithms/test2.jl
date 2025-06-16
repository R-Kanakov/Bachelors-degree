using NonlinearSolve
using DynamicalSystems
using ForwardDiff
using DifferentialEquations
using LinearAlgebra

function newton_section(ds, u0, T, Δt)
    D = length(u0)

    function phase_velocity(u)
        return ds.integ.f.f(u, ds.p0, 0.0)
    end

    v = phase_velocity(u0)

    inds = findall(x -> abs(x) > 1e-8, v)
    if length(inds) < 2
      println(v, u0)
      error("Не удалось найти две ненулевые компоненты фазовой скорости")
    end
    i, j = inds[1], inds[2]

    println("Выбраны компоненты секущей: ", i, ", ", j)

    f = (err, x, p) -> begin
        if isinplace(ds) 
          u = @view x[1:D]
        else
          u = SVector{D}(x[1:D])
        end
        T_local = x[end]

        tspan = (0.0, T_local)
        save_times = 0.0:Δt:T_local

        prob = remake(ds.integ.sol.prob; u0 = u, tspan = tspan)
        sol = solve(prob; saveat = save_times, alg = TRBDF2(), maxiters=1_000_000)

        diff_vec = sol.u[end] .- sol.u[1]

        err[1] = diff_vec[i]
        err[2] = diff_vec[j]

        err[3] = T_local - sol.t[end]
    end

    x0 = [u0; T]

    resid_prototype = zeros(3)

    prob = NonlinearLeastSquaresProblem(
        NonlinearFunction(f, resid_prototype = resid_prototype), x0)

    sol = solve(prob, NonlinearSolve.NewtonRaphson())

    return sol
end


#=
function newton(ds, u0, T, Δt)
  D = dimension(ds)

  f = (err, v, p) -> begin
    if isinplace(ds) 
      u0 = @view v[1:D]
    else
      u0 = SVector{D}(v[1:D])
    end
    T = v[end]

    tspan = (0.0, T)

    T_float = ForwardDiff.value(T)

    save_times = 0.0:Δt:T_float

    sol = solve(SciMLBase.remake(ds.integ.sol.prob; u0 = u0, tspan = tspan);
                DynamicalSystemsBase.DEFAULT_DIFFEQ..., ds.diffeq..., saveat = save_times,
                alg=Rodas5(), maxiters=1_000_000)

    for i in 1:D
        err[i] = sol.u[end][i] - sol.u[1][i]
    end

    err[end] = T - sol.t[end]
  end

  prob = NonlinearLeastSquaresProblem(
      NonlinearFunction(f, resid_prototype = zeros(4)), [u0..., T])

  return solve(prob, NonlinearSolve.NewtonRaphson())
end
=#

function main()
  system = Systems.lorenz()

  u0 = [1.0, 2.0, 3.0]
  T  = 6.0200000000000005
  Δt = 0.01

  sol = newton_section(system, u0, T, Δt)

  println("Решение: ", sol.u)
end

main()
