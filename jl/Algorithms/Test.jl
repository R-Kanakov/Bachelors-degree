using DynamicalSystems
using LinearAlgebra

include("Newton.jl")
include("LM.jl")

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
