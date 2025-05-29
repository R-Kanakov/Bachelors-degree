using GLMakie

function main()
  ε        = 0.4
  length   = 10
  n_points = 100

  x      = range(-length, stop = length, length = n_points)
  y_main = x .^ 2 * 0.005
  y_dash = .-(x .^ 2 * 0.001)
  z      = zeros(n_points)

  l, r = -2, 2

  θ = LinRange(0, 2π, 100)

  filtered_indices = findall(l .<= x .<= r)

  len = size(filtered_indices, 1)

  t = range(l, r, len)

  X = [x[filtered_indices[i]] + t[i] for θ_i in θ, i in 1:len]
  Y = [y_main[filtered_indices[i]] + ε * cos(θ_i) for θ_i in θ, i in 1:len]
  Z = [ε * sin(θ_i) for θ_i in θ, t_i in t]

  fig = Figure(size = (800, 600))
  ax  = Axis3(fig[1, 1])

  lines!(ax, x, y_main, z, color = :black, linewidth = 2)
  lines!(ax, x, y_dash .- ε * 0.85, z, color = :black, linestyle = :dash)

  lines!(ax, x[filtered_indices],
             y_main[filtered_indices],
             z[filtered_indices],
         color = :red, linewidth = 2)

  lines!(ax, x[filtered_indices],
             y_dash[filtered_indices] .- ε * 0.85,
             z[filtered_indices],
         color = :red, linestyle = :dash)

  surface!(ax, X, Y, Z,
           color = fill(:gray, size(X)...),
           alpha = 0.5, shading = NoShading)

  save_path = joinpath(dirname(@__FILE__), "3.png")
  save(save_path, fig)
end

main()
