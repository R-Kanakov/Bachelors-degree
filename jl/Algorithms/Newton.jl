using NonlinearSolve
using DynamicalSystems
using ForwardDiff
using DifferentialEquations
using LinearAlgebra
using StaticArrays

function newton_method(ds::CoupledODEs, u0::U where {U<:AbstractArray{<:Real}}, T::Float64, Δt::Float64, t_min::Float64)
    # Размерность системы
    D = dimension(ds)

    # v - значение вектора фазовой скорости в точке u0
    v = ds.integ.f.f(u0, ds.p0, 0.0)

    # Поиск ненулевых компонент для построения секущей
    inds = findall(x -> abs(x) > 1e-8, v)
    #length(inds) < 2 && (println("Newton: Too many zeros"); return nothing)

    # maximum(inds)

    # Вычисление номеров компонент, в которых располагается секущая
    i, j = inds[1], inds[2]
    # Если i = 1, j = 2, нужно запомнить, что компонента 3 не будет изменяться
    # Обозначим её как k
    k = filter((x) -> (return x != i && x != j), [i for i in 1:D])[1]

    #@assert(length(inds) < 3)
    println("len(inds) = $(length(inds))")

    # Определение невязки
    f = (err, x, p) -> begin
      # Создание вектора начальных условий
      if isinplace(ds)
        u = @view x[[i, j]]
        u[k] = u0[k]
      else
        u = SVector{D}(
          if i == 1
            if j == 2
              (x[i], x[j], u0[k])
            else
              (x[i], u0[k], x[k])
            end
          elseif j == 1
            if i == 2
              (x[j], x[i], u0[k])
            else
              (x[j], u0[k], x[i])
            end
          else
            if i == 2
              (u0[k], x[i], x[j])
            else
              (u0[k], x[j], x[i])
            end
          end
        )
      end

      # Предположительное время периода
      T_local = x[end]

      # Диапазоны, в которых будем интегрировать систему
      tspan = (0.0, T_local)
      save_times = 0.0:Δt:T_local

      # Интегрирование системы для начального условия u на время T_local
      sol = solve(remake(ds.integ.sol.prob; u0 = u, tspan = tspan); 
                  DynamicalSystemsBase.DEFAULT_DIFFEQ..., ds.diffeq..., saveat = save_times)

      # Разница между конечным и начальным состоянием
      diff_vec = sol.u[end] .- sol.u[1]

      # Формирование ошибки
      err[1] = diff_vec[i]
      err[2] = diff_vec[j]
      err[3] = abs(sol.t[end] - sol.t[1])
    end

    x0 = [u0; T]

    prob = NonlinearLeastSquaresProblem(NonlinearFunction(f, resid_prototype = zeros(3)), x0)

    sol = solve(prob, NonlinearSolve.NewtonRaphson())

    u0_new = sol.u[1:end-1]
    T_new  = sol.u[end]

    T_new < t_min && (println("Newton: Period $T_new is less than t_min"); return nothing)

    return (state = trajectory(ds, T_new - Δt, u0_new; Δt = Δt)[1], period = T_new)
end
