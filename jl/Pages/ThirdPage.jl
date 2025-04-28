include("../poincare.jl")
include("../util.jl")


function initThirdPage(screen, fig, u0s, total_time, system, psos, pageClosedT, i_orbits)
  current_size = Makie.size(fig.scene)
  delete_all!(fig)

  trs = [trajectory(system[], total_time[], u0, Δt = 0.01) for u0 ∈ u0s[]]
  trs = filter_orbits(trs)
  for _ in eachindex(u0s[])
    i_orbits[] = push!(i_orbits[], true)
  end

  D = dimension(system[])
  
  psos[] = find_best_poincaresos(trs, D, total_time[])
  
  if D == 4
    if typeof(trs) <: Vector{Tuple{S, T}} where {S<:AbstractStateSpaceSet, T<:StepRangeLen}
        trs = [trs[i][1][:, 1:3] for i in eachindex(trs)]
    else
        trs = [trs[i][:, 1:3] for i in eachindex(trs)]
    end

    psos_ = []

    for k in eachindex(psos[])
      psos_j = []
      for j in 1:3
        push!(psos_j, (s = psos[][k][j].s[:, 1:3], v = psos[][k][j].v))
      end
      push!(psos_, psos_j)
    end
    
    psos[] = psos_
  end

  fig, psos = scan_poincaresos(trs, screen, current_size, psos, i_orbits, total_time[], pageClosedT; linekw = (transparency = true,))
end
