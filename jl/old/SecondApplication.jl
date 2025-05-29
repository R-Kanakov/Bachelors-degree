using CairoMakie
using LaTeXStrings

function main()
  fig = Figure(size = (600, 200))

  ax1 = Axis(fig[1, 1], title = L"l_1")
  ax2 = Axis(fig[1, 2], title = L"l_2")
  axi = Axis(fig[1, 3], title = L"l_\infty")

  for ax in [ax1, ax2, axi]
    limits!(ax, -2, 2, -2, 2)
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
  end

  θ = LinRange(0, 2π, 1000)

  l1_x = [sign(cos(t)) * abs(cos(t))^2 for t in θ]
  l1_y = [sign(sin(t)) * abs(sin(t))^2 for t in θ]
  lines!(ax1, l1_x, l1_y)
  scatter!(ax1, 0, 0, markersize = 5)

  l2_x = [cos(t) for t in θ]
  l2_y = [sin(t) for t in θ]
  lines!(ax2, l2_x, l2_y)
  scatter!(ax2, 0, 0, markersize = 5)

  li_x = [1, 1, -1, -1, 1]
  li_y = [1, -1, -1, 1, 1]
  lines!(axi, li_x, li_y)
  scatter!(axi, 0, 0, markersize = 5)

  save_path = joinpath(dirname(@__FILE__), "2.png")
  save(save_path, fig)
end

main()