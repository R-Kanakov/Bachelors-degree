using DynamicalSystems
using LinearAlgebra

include("Newton.jl")
include("LM.jl")

function main()
  system = Systems.chua()

  ig = [1.0883306197960876, -0.03702073227722233, -0.22105768451157518]
  T  = 6.48
  Δt = 0.01

  res = trajectory(system, T - Δt, ig; Δt = Δt)[1]
  Δ   = LinearAlgebra.norm(res[end] .- res[1])

  res_newton = newton_method(system, ig, T, Δt, 1.0)
  !isnothing(res_newton) && (Δ_newton = LinearAlgebra.norm(res_newton.state[end] .- res_newton.state[1]))

  res_lm = lm(system, ig, T, Δt, 1.0)
  !isnothing(res_lm) && (Δ_lm = LinearAlgebra.norm(res_lm.state[end] .- res_lm.state[1]))

  println("Δ = $Δ, T = $T")
  !isnothing(res_newton) && (println("Newton begin state = $(res_newton.state[1])"))
  !isnothing(res_newton) && (println("Δ_newton = $Δ_newton, T = $(res_newton.period)"))
  !isnothing(res_lm) && (println("LM begin state = $(res_lm.state[1])"))
  !isnothing(res_lm) && (println("Δ_lm = $Δ_lm, T = $(res_lm.period)"))
end

main()
