using DynamicalSystems
using CairoMakie
using LaTeXStrings

function main()
  Lorenz = Systems.lorenz()

  u0 = [10, 10, 10.9]

  total_time = 15

  range = -0.001:0.0005:0.001

  tr = [trajectory(Lorenz, total_time, u0 .+ ϵ; Dt = 0.01) for ϵ in range]

  CairoMakie.activate!()
  fig = Figure(size = (919, 237))
  ax  = Axis(fig[1, 1]; xlabel = L"\mathbf{t}", ylabel = L"\mathbf{z}")

  for i in 1:length(range)
    lines!(ax, tr[i][2], columns(tr[i][1])[3])
  end

  save_path = joinpath(dirname(@__FILE__), "1.png")
  save(save_path, fig)
end

main()