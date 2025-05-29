using DynamicalSystems, OrdinaryDiffEq
using GLMakie
using LaTeXStrings

include("../util.jl")


function create_u_boxes(grid, i, u0s, D)
  latex_expr = latexstring("u_{", i, "} = ")
  Label(grid[i, 1], latex_expr)

  for j in 1:D
    tb = Textbox(grid[i, j + 1], placeholder = string(u0s[][i][j]), validator = Float64)
    on(tb.stored_string) do u u0s[][i][j] = parse(Float64, u) end
  end
end

function add_u0(u0s, grids, D)
  len = length(u0s[])
  if len < 6
    insertrows!(grids[][3], len + 1, 1)
    u0s[] = push!(u0s[], rand(Float64, D) .% 1)
    create_u_boxes(grids[][3], len + 1, u0s, D)
  end
end

function del_u0(u0s, grids, D)
  len = length(u0s[])

  content = contents(grids[][3])
  filter!(elt->(return !(typeof(elt) <: Button)), content)
  len_contents = length(content)

  if len > 1
    for i in 0:D
      delete!(content[len_contents - D + i])
    end
  end
  trim!(grids[][3])
  
  u0s[] = u0s[][1:len - 1]
end

function create_button_next(fig, grid, flag)
  buttonNext = Button(grid[1, 3], label = "Next",
                      halign = :right, valign = :top,
                      tellwidth = false, tellheight = false)
  on(buttonNext.clicks) do _
    delete_all!(fig)
    flag[] = true
  end
end


function lorenz_u0s(u0s)
  u0s[] = push!(u0s[], [1.0, 1.0, 1.0])
  u0s[] = push!(u0s[], [-1.0, -1.31596, -1.354461])
  u0s[] = push!(u0s[], [0.0, -0.25, 0.42081])
end

function henonheiles_u0s(u0s)
  u0s[] = push!(u0s[], [0.0, -0.25, 0.42081, 0.0])
  u0s[] = push!(u0s[], [0.0, 0.1, 0.5, 0.0])
  u0s[] = push!(u0s[], [0.0, -0.31596, 0.354461, 0.0591255])
end

function rikitake_u0s(u0s)
  u0s[] = push!(u0s[], [1.0, 1.0, 1.0])
  u0s[] = push!(u0s[], [-1.0, -1.31596, -1.354461])
  u0s[] = push!(u0s[], [0.0, -0.25, 0.42081])
end

function gissinger_u0s(u0s)
  u0s[] = push!(u0s[], [1.0, 1.0, 1.0])
  u0s[] = push!(u0s[], [-1.0, -1.31596, -1.354461])
  u0s[] = push!(u0s[], [0.0, -0.25, 0.42081])
end

function roessler_u0s(u0s)
  u0s[] = push!(u0s[], [1.0, 1.0, 1.0])
  u0s[] = push!(u0s[], [-1.0, -1.31596, -1.354461])
  u0s[] = push!(u0s[], [0.0, -0.25, 0.42081])
end

function chua_u0s(u0s)
  u0s[] = push!(u0s[], [0.0, -0.25, 0.42081])
  u0s[] = push!(u0s[], [1.0, 1.0, 1.0])
  u0s[] = push!(u0s[], [-1.0, -1.31596, -1.354461])
end


function initSecondPage(fig, ds, rule, flag, u0s, total_time, system)
  grids = Observable([GridLayout(fig[i, 1], tellwidth = false, tellheight = false) for i in 1:3])

  Label(grids[][1][1, 2], "Enter initial states", 
                    halign = :center, valign = :top,
                    tellwidth = false, tellheight = false)
  create_button_next(fig, grids[][1], flag)
  rowsize!(grids[][1], 1, Relative(1/5))

  D = 0

  if ds != "Enter"
    ds[] == "Lorenz"      && (system[] = Systems.lorenz();      lorenz_u0s(u0s))
    ds[] == "Henonheiles" && (system[] = Systems.henonheiles(); henonheiles_u0s(u0s))
    ds[] == "Rikitake"    && (system[] = Systems.rikitake();    rikitake_u0s(u0s))
    ds[] == "Gissinger"   && (system[] = Systems.gissinger();   gissinger_u0s(u0s))
    ds[] == "Roessler"    && (system[] = Systems.roessler();    roessler_u0s(u0s))
    ds[] == "Chua"        && (system[] = Systems.chua();        chua_u0s(u0s))
    D = dimension(system[])
  else
    # Definition of DS via parsed rule...
    system[] = Systems.lorenz()
    D = dimension(system[])

    if D == 3
      lorenz_u0s(u0s)
    elseif D == 4
      henonheiles_u0s(u0s)
    end
  end

  Label(grids[][2][1, 1], "Evaluation time = ")

  tb = Textbox(grids[][2][1, 2], placeholder = string(total_time[]), validator = Int64)
  on(tb.stored_string) do t total_time[] = parse(Int, t) end

  rowsize!(grids[][2], 1, Relative(1/5))

  for i in 1:3 create_u_boxes(grids[][3], i, u0s, D) end

  butCreate = Button(grids[][3][4, 2], label = "+")
  on(butCreate.clicks) do _ add_u0(u0s, grids, D) end

  butCreate = Button(grids[][3][4, 3], label = "-")
  on(butCreate.clicks) do _ del_u0(u0s, grids, D) end
end
