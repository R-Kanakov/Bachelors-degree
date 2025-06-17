using LaTeXStrings
using CairoMakie, GLMakie
using SparseArrays
using Distances: evaluate
using LinearAlgebra
using Base.Filesystem
using LightGraphs

include("../Algorithms/LM.jl")
include("../Algorithms/Newton.jl")

function callback(system, psos, v, j, dir, u0s, t_max, t_min, ε, δ, Δt, total_time, saveR, saveO, ds)
  R, res = [], []

  current_dir = dirname(dirname(dirname(@__FILE__)))

  logs_path = joinpath(current_dir, "logs", "$ds")

  !isdir(logs_path) && mkdir(logs_path)

  file_prefix = "$(ds)_"
  file_extension = ".txt"

  new_file_path = create_next_file(logs_path, file_prefix, file_extension)

  log = open(new_file_path, "w")

  # Logging initial state
  write(log, "t_max = $(t_max[]), t_min = $(t_min[])\n")
  write(log, "ε = $(ε[]), δ = $(δ[])\n\n")
  write(log, "Initial states are\n")
  for i in eachindex(u0s)
    write(log, "    u[$i] = $(u0s[i])\n")
  end
  write(log, "\nNumber of founded dots is $(sum(length(psos_) for psos_ in psos))\n")
  write(log, "\nPoincare surface of sections defined by:\n")
  write(log, " - Direction of movement along the axis $(('x':'z')[j])\n")
  write(log, " - Value in that direction is $v\n")
  write(log, " - Direction of plane intersection comes with ") 
  write(log, "$(if dir == true "increasing coordinate (towards the positive normal)" 
                else "decreasing coordinates (towards the negative normal)"
                end)")
  write(log, "")
  write(log, "\n\n")

  # Making paths to Png and Mkv directories
  if saveR[]
    png_path = joinpath(current_dir, "Png", "$ds")
    rm(png_path, recursive = true, force = true)
    mkdir(png_path)
  end

  if saveO[]
    mkv_path = joinpath(current_dir, "Mkv", "$ds")
    rm(mkv_path, recursive = true, force = true)
    mkdir(mkv_path)
  end

  # Integrating each dot on psos
  # psos - вектор, элементом которого является набор точек
  # из сечения Пуанкаре для определённого начального условия

  for i in eachindex(psos) # идём по каждому начальному значению
    push!(R, [])

    # Интегрируем систему на время t_max для каждой точки из плоскости Пуанкаре
    trs = [trajectory(system, t_max[], u0, Δt = Δt) for u0 in psos[i]]

    # Для каждой траектории смотрим в первый столбец, если существует возвращаение
    # в ε окрестность, то сохраняем период
    Rs = [let d = distance_col_s(tr[1], ε[], Int(t_min[] / Δt))
            if d.b == true
              r = d.r
            end
          end for tr in trs]
    
    for k in eachindex(Rs)
      if !isnothing(Rs[k])
        push!(R[i], [trs[k][1][1], Rs[k]]) # R[начальное условие][все точки для которых есть возвращение]
      end
    end
  end

  (all(length(R[i]) == 0 for i in eachindex(psos))) &&
    (write(log, "0 orbits found.\nTry to increase total_time, t_max or ε.");
     close(log);
     return)

  lm_u0s, lm_t0s = [], []

  for i in eachindex(R) # for each u0s
    for j in eachindex(R[i]) # for each dot in u0s
      push!(res, lm(system, collect(R[i][j][1]), R[i][j][2] * Δt, Δt, t_min[]))
      push!(lm_u0s, R[i][j][1])
      push!(lm_t0s, R[i][j][2] * Δt)
    end
  end

  all(length(res[i]) == 0 for i in eachindex(res)) &&
    (write(log, "0 orbits found\nTry to increase total_time, t_max or ε.");
     close(log);
     return)

  psos_res = []

  for i in eachindex(res)
    if !isnothing(res[i])
      write(log, "\nOrbit $i\n")

      Δ = norm(res[i].state[end] .- res[i].state[1])

      write(log, "Before Levenberg-Marquardt:\n",
                 "  Begin state = $(lm_u0s[i])\n",
                 "  T = $(lm_t0s[i])\n",
                 "After Levenberg-Marquardt:\n",
                 "  Begin state = $(res[i].state[1])\n",
                 "  End state   = $(res[i].state[end])\n",
                 "  Δ = $Δ\n",
                 "  T = $(res[i].period)\n")

      Δ > ε[] && (write(log, "Distance has been increased. Skipping orbit\n"); continue)

      try
        psos_ = poincaresos(res[i].state, (j, v); direction = -1, warning = false)
        push!(psos_res, (psos = psos_, i = i))
      catch err
        write(log, "Psos is empty\n")
        continue
      end

      write(log, "Psos intersection number is $(length(psos_res[end].psos))\n")
    end
  end

  write(log, "\nOrbits left after LM/Newton is $(length(psos_res))\n")

  psos_res = group_psos(psos_res)

  d = 3
  if j == 1
    d = 2
  elseif j == 2 || j == 3
    d = 1
  end

  animate = []

  for i in eachindex(psos_res)
    len_i = length(psos_res[i])

    write(log, "\n$(length(psos_res[i][1].psos)) intersections have orbits with indexes\n")
    write(log, "[i = $(psos_res[i][1].i)")

    for j in 2:len_i
      write(log, "\n i = $(psos_res[i][j].i)")
    end

    write(log, "]\n")

    if len_i != 1 # if group contain several orbits
      psos = []

      for j in eachindex(psos_res[i])         # For each psos in group
        k = min_index(psos_res[i][j].psos, d) # We find minimal component
        if k != 1                             # And if it not first
                                              # We rotate psos so it will be first
          push!(psos, vcat(psos_res[i][j].psos[k:end], 
                      psos_res[i][j].psos[1:k - 1]))
        elseif k == 1                         # Else everything already in its place
          push!(psos, psos_res[i][j].psos)
        end
      end

      # Grouping our psoses by distance δ
      grouped_ind = compare_psoses(psos, δ[])

      if length(grouped_ind) == 1
        write(log, "All orbits were same, we'll consider first of them\n")
        saveO[] && (append!(animate, psos_res[i][1].i))
      elseif length(grouped_ind) != len_i
        write(log, "Orbits with indexes\n  ")

        filtered_group     = filter(x -> (length(x) != 1), grouped_ind)
        filtered_group_len = length(filtered_group)

        for j in 1:filtered_group_len
          if j != filtered_group_len
            write(log, "$(map(x -> x.i, psos_res[i][[filtered_group[j]...]])),\n  ")
          else
            write(log, "$(map(x -> x.i, psos_res[i][[filtered_group[j]...]]))\n")
          end
        end

        write(log, "were same\n")

        if saveO[] || saveR[]
          for j in eachindex(grouped_ind)
            append!(animate, psos_res[i][grouped_ind[j][1]].i)
          end
        end
      else
        write(log, "All orbits were different\n")
        if saveO[] || saveR[]
          for j in 1:len_i
            append!(animate, psos_res[i][j].i)
          end
        end
      end
    else # if group contain single orbit
      if saveO[] || saveR[]
        append!(animate, psos_res[i][1].i)
      end
    end
  end

  if saveO[]
    GLMakie.activate!()
    for i in animate
      min, max = minmaxima(res[i].state)

      save_path = joinpath(mkv_path, "Orbit_$(i).mp4")

      len  = length(res[i].state)
      time = Observable(0.0)
      k    = @lift(round(Int, clamp(len * $time, 1, len)))
      azim = Observable(1.275)
      z    = @lift(clamp(1.0 * $azim, 1.275, 2pi + 1.275))

      figi = Figure(size = (1000, 1000))

      axi  = Axis3(figi[1, 1], azimuth = @lift($z * pi))
      limits!(axi, min[1], max[1],
                   min[2], max[2],
                   min[3], max[3])

      dynamic_mas = @lift(res[i].state[1:$k])
      lines!(axi, dynamic_mas, color = RGBA{Float32}(0.0, 0.44705883, 0.69803923, 1.0), linewidth = 3)

      eval_range    = range(0, 1, step = 1 / len)
      rotate_range  = range(1, 3, step = 10e-3)
      rotate_range_ = range(3, 5, step = 10e-3)
      timestamps    = vcat(eval_range, rotate_range, rotate_range_) 

      record(figi, save_path, timestamps) do t
        if t <= 1
          time[] = t
        elseif t < 3
          azim[] = t - 1 + 1.275
        elseif t == 3
          i1, i2 = find_element_with_i(psos_res, i)
          scatter!(axi, psos_res[i1][i2].psos, color = :red, markersize = 10)

          ss = collect(copy(min))
          ws = collect(copy(max) .- copy(min))

          ws[j] = 0
          ss[j] = v

          o = Makie.Point3f(ss...)
          w = Makie.Point3f(ws...)

          mesh!(axi, Makie.Rect3f(o, w); color = RGBAf(0.2, 0.2, 0.25, 0.5), transparency = true)
        elseif t <= 5
          azim[] = t - 1 + 1.275
        end
      end
    end
  end

  if saveR[]
    GLMakie.activate!()
    Y = [trajectory(system, total_time, u0; Dt = Δt)[1] for u0 in u0s]

    min, max, _ = total_minmaxima(Y)

    for i in animate
      save_path = joinpath(png_path, "Orbit_$(i).png")
      figi = Figure(size = (2000, 1000))
      axi1 = Axis3(figi[1, 1], azimuth = 0.5 * pi, titlesize = 40, title = "Original attractor")
      limits!(axi1, min[1], max[1],
                    min[2], max[2],
                    min[3], max[3])
      for j in eachindex(Y)
        lines!(axi1, Y[j][:, 1], Y[j][:, 2], Y[j][:, 3]; linewidth = 1, color = RGBA{Float32}(0.0, 0.44705883, 0.69803923, 1.0))
      end

      axi2 = Axis3(figi[1, 2], azimuth = 0.5 * pi, titlesize = 40, title = "Reconstructed attractor")
      limits!(axi2, min[1], max[1],
                    min[2], max[2],
                    min[3], max[3])
      lines!(axi2, res[i].state[:, 1], res[i].state[:, 2], res[i].state[:, 3], color = :blue)

      save(save_path, figi)

      save_path = joinpath(png_path, "Orbit_$(i)_1.png")
      figi = Figure(size = (1000, 1000))
      axi1 = Axis3(figi[1, 1], azimuth = 0.5 * pi, titlesize = 40, title = "Original attractor with periodic orbit")
      limits!(axi1, min[1], max[1],
                    min[2], max[2],
                    min[3], max[3])

      for j in eachindex(Y)
        lines!(axi1, Y[j][:, 1], Y[j][:, 2], Y[j][:, 3]; linewidth = 1, color = RGBA{Float32}(0.0, 0.44705883, 0.69803923, 0.4))
      end

      lines!(axi1, res[i].state[:, 1], res[i].state[:, 2], res[i].state[:, 3]; linewidth = 8, color = :red)

      save(save_path, figi)
    end
  end

  close(log)
end

function create_button_next(fig, callback, args...)
  buttonNext = Button(fig[1, 2], label = "Next")
  on(buttonNext.clicks) do _
    call = @async callback(args...)
    delete_all!(fig)

    try
      wait(call)
    catch error
      showerror(stderr, error)
      println(stderr)
      @error "An error occurred" exception=(error, catch_backtrace())
    end
  end
end

function distance_col_s(x::SSSet, ε::Real, t_min::Int)
  for i in t_min:length(x)
      @inbounds if evaluate(Euclidean(), x[1], x[i]) ≤ ε
        return (b = true, r = i)
      end
  end
  return (b = false, )
end

function find_element_with_i(psos_res, I)
  for (outer_index, inner_array) in enumerate(psos_res)
      for (inner_index, element) in enumerate(inner_array)
          if element.i == I
              return (outer_index, inner_index)
          end
      end
  end
  return nothing
end

function group_by_key(arr, keyfunc)
  groups = Dict{Any, Vector{Any}}()
  for el in arr
      k = keyfunc(el)
      push!(get!(groups, k, Vector{Any}()), el)
  end
  return collect(values(groups))
end

function group_psos(psos_res)
  return group_by_key(psos_res, x -> length(x.psos))
end

function min_index(data::SSSet, dir)
  mi = Vector(data[1])

  res = [1, 1, 1]

  for i in eachindex(data)
    if data[i][dir] < mi[dir]
        mi[dir]  = data[i][dir]
        res[dir] = i
    end
  end

  return res[dir]
end

function compare_psoses(psoses, δ)
  differ = []

  len_psos = length(psoses[1])
  n        = length(psoses)

  for i in 1:len_psos
    for j in 1:n
      for k in j + 1:n
        if evaluate(Euclidean(), psoses[j][i], psoses[k][i]) >= δ
          push!(differ, (j, k))
        end
      end
    end
  end

  unique!(differ)

  return find_optimal_groups(n, differ)
end

function find_optimal_groups(n, differ)
  forbidden = Set((min(u,v), max(u,v)) for (u,v) in differ)

  g = SimpleGraph(n)
  for u in 1:n
      for v in u+1:n
          if !( (u,v) ∈ forbidden )
              add_edge!(g, u, v)
          end
      end
  end

  components = connected_components(g)

  result = []
  for comp in components
      group = tuple(comp[1], comp[2:end]...)
      push!(result, group)
  end

  return result
end

function create_next_file(dir_path, file_prefix, file_extension)
  files = readdir(dir_path)

  filtered_files = filter(file -> startswith(file, file_prefix) && endswith(file, file_extension), files)

  numbers = Int[]
  for file in filtered_files
      suffix = file[length(file_prefix)+1:end-length(file_extension)]
      try
          num = parse(Int, suffix)
          push!(numbers, num)
      catch end
  end

  max_number = isempty(numbers) ? 0 : maximum(numbers)

  new_file_number = max_number + 1
  new_file_name = "$(file_prefix)$(new_file_number)$(file_extension)"
  new_file_path = joinpath(dir_path, new_file_name)

  touch(new_file_path)

  return new_file_path
end

function to_line(tuple)
  field_names = fieldnames(typeof(tuple))

  line_parts = String[]
  for name in field_names
      value = getproperty(tuple, name)
      push!(line_parts, "$(name) = $(value)")
  end

  line = "(" * join(line_parts, ", ") * ")"
  return line
end


function initFourthPage(screen, fig, system, final_psos, i_orbits, ds, u0s, total_time)
  current_size = Makie.size(fig.scene)
  figure = Figure(size = current_size)
  display(screen, figure)

  psos = [if i_orbits[][i] == true 
            final_psos[].psos[i]
          else
            nothing
          end for i in eachindex(i_orbits[])]

  u0s = [if i_orbits[][i] == true 
            u0s[i]
          else
            nothing
          end for i in eachindex(i_orbits[])]
  filter!(!isnothing, psos)
  filter!(!isnothing, u0s)

  v    = final_psos[].v
  j    = final_psos[].j
  dir  = final_psos[].d

  saveR = Observable(false)
  saveO = Observable(false)
  t_max = Observable(30.0)
  t_min = Observable(1.0)
  ε     = Observable(0.01)
  δ     = Observable(0.05)
  Δt    = 0.01

  ds == "Chua" && (ε[] = 0.005)

  label  = Label(figure[1, 1], L"\text{Enter}\,\epsilon,\,\text{maximum and minimum period}")

  label1 = Label(figure[2, 1], L"t_{max}=")
  tb1    = Textbox(figure[2, 2], placeholder = string(t_max[]), validator = Float64)

  # TODO: make denominator for tb2 (t_min > 0)
  label2 = Label(figure[3, 1], L"t_{min}=")
  tb2    = Textbox(figure[3, 2], placeholder = string(t_min[]), validator = Float64)

  # TODO: make denominator for tb3 (ε > 0)
  label3 = Label(figure[4, 1], L"\epsilon=")
  tb3    = Textbox(figure[4, 2], placeholder = string(ε[]),     validator = Float64)

  # TODO: make denominator for tb4 (δ > 0)
  label4 = Label(figure[5, 1], L"\delta=")
  tb4    = Textbox(figure[5, 2], placeholder = string(δ[]),     validator = Float64)

  label5 = Label(figure[6, 1], "Save orbits picture? ")
  cb1    = Checkbox(figure[6, 2], checked = false)

  label6 = Label(figure[7, 1], "Save orbits animation? ")
  cb2    = Checkbox(figure[7, 2], checked = false)

  on(tb1.stored_string) do t t_max[] = parse(Float64, t) end
  on(tb2.stored_string) do t t_min[] = parse(Float64, t) end
  on(tb3.stored_string) do e ε[]     = parse(Float64, e) end
  on(tb4.stored_string) do d δ[]     = parse(Float64, d) end
  on(cb1.checked)       do s saveR[] = s                 end
  on(cb2.checked)       do s saveO[] = s                 end

  create_button_next(figure, callback, system[], psos, v, j, dir, u0s, t_max,
                     t_min, ε, δ, Δt, total_time, saveR, saveO, ds)
end
