using LaTeXStrings
using CairoMakie, GLMakie
using SparseArrays
using Distances: evaluate
using LinearAlgebra

include("../Algorithms/LM.jl")
include("../util.jl")


# TODO:
# 1. Correct log writing
# 2. Make size and markersize depent from s (t_max)
# 3. For now we watch only first return
# 4. Correct orbits selection
# 5. If D == 4 algorithm doesn't work. There has to be a choice as to which component we project onto


function callback(system, psos, t_max, t_min, ε, Δt, saveR, saveO, ds)
  R, res, indices = [], [], []

  if saveR[]
    current_dir = dirname(dirname(dirname(@__FILE__)))
    png_path = joinpath(current_dir, "Png", "$ds")
    rm(png_path, recursive = true, force = true)
    mkdir(png_path)
  end

  if saveO[]
    current_dir = dirname(dirname(dirname(@__FILE__)))
    mkv_path = joinpath(current_dir, "Mkv", "$ds")
    rm(mkv_path, recursive = true, force = true)
    mkdir(mkv_path)
  end

  for i in eachindex(psos[])
    push!(R, [])
    for j in 1:3
      trs = [trajectory(system, t_max[], u0, Δt = Δt) for u0 in psos[][i][j].s]

      if saveR[]
        Rs = [let d = distance_col_s(tr[1], ε[], Int(t_min[] / Δt))
                if d.b == true
                  (m = sparse(RecurrenceAnalysis.tril(RecurrenceMatrix(tr[1], ε[]).data, -1)), r = d.r)
                end
              end for tr in trs]
      else
        Rs = [let d = distance_col_s(tr[1], ε[], Int(t_min[] / Δt))
                if d.b == true
                  (r = d.r, )
                end
              end for tr in trs]
      end

      if saveR[]
        for k in eachindex(Rs)
          if !isnothing(Rs[k])
            x, y = coordinates(Rs[k].m)
            s = size(Rs[k].m, 1)

            CairoMakie.activate!()

            figR = Figure(size = (5000, 5000))
            axR  = Axis(figR[1, 1], titlesize = 50, title = "Recurrence matrix")
            CairoMakie.scatter!(axR, x, y; color = :red, markersize = 5)
            CairoMakie.lines!(axR, [0, s], [0, s]; color = :black)
            axR.limits = ((1, s), (1, s))
            axR.aspect = 1

            save_path = joinpath(png_path, "Recurrence_Matrix_$(i)-Orbit_$(j)-PsosDir_$(k)-Dot.png")
            save(save_path, figR)
          end
        end
      end

      for k in eachindex(Rs)
        if !isnothing(Rs[k])
          push!(R[i], [trs[k][1][1], j, Rs[k].r])
          push!(indices, (j = j, v = psos[][i][j].v))
        end
      end
    end
  end

  (all(length(R[i]) == 0 for i in eachindex(psos[]))) && (@error "0 orbits found\nTry to increase total_time, t_max or ε.")

  for i in eachindex(R)
    for j in eachindex(R[i])
      push!(res, let lm = lm(system, collect(R[i][j][1]), R[i][j][3] * Δt, Δt, t_min[])
                   !isnothing(lm) && (lm)
                 end
           )
    end
  end

  filter!(elt->(return !(elt == false)) , res)

  (all(length(res[i]) == 0 for i in eachindex(res))) && (@error "0 orbits found\nTry to increase total_time, t_max or ε.")

  psos_res = []

  for i in eachindex(res)
    println("\nOrbit $i")

    Δ = norm(res[i].state[end] .- res[i].state[1])

    println("Begin state = ", res[i].state[1], "\nEnd state   = ", res[i].state[end], "\nΔ = ", Δ, "\nT = ", res[i].period)

    println("Psos direction and value = ", indices[i])

    Δ > ε[] && (println("Distance has been increased. Skipping orbit"); continue)

    try
      psos_ = poincaresos(res[i].state, (indices[i].j, indices[i].v); direction = -1, warning = false)
      push!(psos_res, (psos = psos_, i = i, direction = indices[i].j, v = indices[i].v))
    catch err
      println("Psos = empty")
      continue
    end

    println("Psos intersection number is ", length(psos_res[end].psos))

    if saveO[]
      min, max = minmaxima(res[i].state)

      save_path = joinpath(mkv_path, "Orbit$i.mp4")

      len  = length(res[i].state)
      time = Observable(0.0)
      k    = @lift(round(Int, clamp(len * $time, 1, len)))

      GLMakie.activate!()

      figi = Figure(size = (1000, 1000))
      axi  = Axis3(figi[1, 1])
      limits!(axi, min[1], max[1], min[2], max[2], min[3], max[3])

      dynamic_mas = @lift(res[i].state[1:$k])
      lines!(axi, dynamic_mas, color = :blue, linewidth = 3)

      timestamps = range(0, 1, step = 1 / len)

      record(figi, save_path, timestamps) do t
        time[] = t
      end
    end
  end

  psos_res = group_psos(psos_res)

  for i in eachindex(psos_res)
    println("\n", psos_res[i])
    if length(psos_res[i]) != 1
      m = graph_construction(psos_res[i])
      println("Distance Matrix:")
      for row in eachrow(m)
        println(row)
      end
    end
  end
end

function create_button_next(fig, callback, args...)
  buttonNext = Button(fig[1, 2], label = "Next")
  on(buttonNext.clicks) do _
    call = @async callback(args...)#
    delete_all!(fig)

    try
      wait(call)
    catch error
      println(error)
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

function group_by_key(arr, keyfunc)
  groups = Dict{Any, Vector{Any}}()
  for el in arr
      k = keyfunc(el)
      push!(get!(groups, k, Vector{Any}()), el)
  end
  return collect(values(groups))
end

function group_psos(psos_res)
  groups_len = group_by_key(psos_res, x -> length(x.psos))
  
  result = Vector{Vector{Any}}()
  
  for group_len in groups_len
      groups_dir = group_by_key(group_len, x -> x.direction)
      
      for group_dir in groups_dir
          groups_v = group_by_key(group_dir, x -> x.v)
          
          append!(result, groups_v)
      end
  end
  
  return result
end

function graph_construction(psoses)
  len = length(psoses)
  m = Matrix{Float64}(undef, len, len)
  m = zero(m)
  for i in eachindex(psoses)
    for j in i + 1:length(psoses)
      distance = set_distance(psoses[i].psos, psoses[j].psos, StrictlyMinimumDistance())
      m[i, j] = m[j, i] = distance
    end
  end
  return m
end

function initFourthPage(screen, fig, system, psos, i_orbits, ds)
  current_size = Makie.size(fig.scene)
  figure = Figure(size = current_size)
  display(screen, figure)

  psos[] = psos[][i_orbits[]]

  saveR = Observable(false)
  saveO = Observable(false)
  t_max = Observable(50.0)
  t_min = Observable(0.5)
  ε     = Observable(0.01)
  Δt    = 0.01

  label  = Label(figure[1, 1], L"\text{Enter}\,\epsilon,\,\text{maximum and minimum period}")

  label1 = Label(figure[2, 1], L"t_{max}=")
  tb1    = Textbox(figure[2, 2], placeholder = string(t_max[]), validator = Float64)

  # TODO: make denominator for tb2 (0.01 <= t_min)
  label2 = Label(figure[3, 1], L"t_{min}=")
  tb2    = Textbox(figure[3, 2], placeholder = string(t_min[]), validator = Float64)
 
  # TODO: make denominator for tb3 (0.05 <= ε <= 1)
  label3 = Label(figure[4, 1], L"\epsilon=")
  tb3    = Textbox(figure[4, 2], placeholder = string(ε[]),     validator = Float64)

  label4 = Label(figure[5, 1], "Save recurrence matrices? ")
  cb1    = Checkbox(figure[5, 2], checked = false)

  label5 = Label(figure[6, 1], "Save orbits animation? ")
  cb2    = Checkbox(figure[6, 2], checked = false)

  on(tb1.stored_string) do t t_max[] = parse(Float64, t) end
  on(tb2.stored_string) do t t_min[] = parse(Float64, t) end
  on(tb3.stored_string) do e ε[]     = parse(Float64, e) end
  on(cb1.checked)       do s saveR[] = s                 end
  on(cb2.checked)       do s saveO[] = s                 end
  
  create_button_next(figure, callback, system[], psos, t_max, t_min, ε, Δt, saveR, saveO, ds)
end
