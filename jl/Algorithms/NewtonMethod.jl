# TODO (Doesn't work right now)

using LinearAlgebra
using DynamicalSystems

function integrate_to_time(system, T)
  t_target = current_time(system) + T
  while current_time(system) < t_target
      step!(system)
      successful_step(system) || break
  end
  return system.integ.u
end

function newton_sec(x1, x2, T, system; tol=1e-8, max_iter=100)
  v = system.integ.f.f(x1, system.p0, 0)
  n = length(x1)
  
  function extended(x)
      reinit!(system, x)
      F_val = integrate_to_time(system, T) - x
      return vcat(F_val, dot(F_val, v))
  end
  
  F1 = extended(x1)
  
  for steps in 1:max_iter
      F2 = extended(x2)
      Δx = x2 - x1
      ΔF = F2 - F1
      
      norm(Δx) < 1e-12 && return (state=x1, period=T)
      
      Δx_conv = convert(Vector{Float64}, Δx)
      ΔF_conv = convert(Vector{Float64}, ΔF)
      
      J = (ΔF_conv * Δx_conv') / dot(Δx_conv, Δx_conv)
      J_pinv = pinv(J)
      
      dx = J_pinv * ΔF_conv
      dx_orth = dx - dot(dx, v) * v / dot(v, v)
      
      x2 .-= dx_orth
      
      if norm(dx_orth) < tol
          println("Converged in $steps steps")
          return (state=x2, period=T)
      end
      
      x1, F1 = copy(x2), copy(F2)
  end
  
  println("Max iterations reached")
  return (state=x2, period=T)
end
