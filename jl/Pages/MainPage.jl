using GLMakie
using DynamicalSystems
using GridLayoutBase

include("../util.jl")


# TODO: Make function that parse strings to DS rule
function parse_strings()
  return true
end

function create_del_button(str, grid)
  delete!(contents(grid[][3])[4])

  str[] = push!(str[], " ")

  tb4 = Textbox(grid[][3][4, 1], placeholder = "Enter 4")
  on(tb4.stored_string) do s str[][4] = s end

  del_btn = Button(grid[][3][5, 1], label = "-")
  on(del_btn.clicks) do _
    delete!(del_btn)
    delete!(tb4)
    trim!(grid[][3])

    str[] = str[][1:3]

    add_btn = Button(grid[][3][4, 1], label = "+")
    on(add_btn.clicks) do _ create_del_button(str, grid) end
  end
end

function create_button_next_main(fig, grids, ds, flag)
  buttonNext = Button(grids[1][1, 3], label = "Next",
                      halign = :right, valign = :top,
                      tellwidth = false, tellheight = false)
  on(buttonNext.clicks) do _
    if ds != "Enter"
      delete_all!(fig)
      flag[] = true
    else
      res = parse_strings()

      if res # Parse Ok
        delete_all!(fig)
        flag[] = true
      else
        @error "Parse error"
      end
    end
  end
end

function initMainPage()
  screen = GLMakie.Screen()
  fig = Figure()

  grids = Observable([GridLayout(fig[i, 1], tellwidth = false, tellheight = false) for i in 1:3])

  Label(grids[][1][1, 2], "Choose Dynamical System", 
                    halign = :center,
                    valign = :top,
                    tellwidth  = false,
                    tellheight = false)
  rowsize!(grids[][1], 1, Relative(1/5))

  menu = Menu(grids[][2][1, 1], 
              options = ["Lorenz", "Henonheiles", "Rikitake", "Gissinger", "Roessler", "Chua", "Enter"],
              default = "Lorenz",
              width   = 200,
              halign  = :center,
              valign  = :top,
              tellwidth  = false,
              tellheight = false)
  rowsize!(grids[][2], 1, Relative(1/5))
  
  ds           = Observable("Lorenz")
  str          = Observable(String[" ", " ", " "])
  flag         = Observable(false)
  created_flag = false
  rule         = nothing

  on(menu.selection) do s
      if s == "Enter"
        if !created_flag
          tb1 = Textbox(grids[][3][1, 1], placeholder = "Enter 1")
          on(tb1.stored_string) do s str[][1] = s end

          tb2 = Textbox(grids[][3][2, 1], placeholder = "Enter 2")
          on(tb2.stored_string) do s str[][2] = s end

          tb3 = Textbox(grids[][3][3, 1], placeholder = "Enter 3")
          on(tb3.stored_string) do s str[][3] = s end

          butCreate = Button(grids[][3][4, 1], label = "+")
          on(butCreate.clicks) do _ create_del_button(str, grids) end
        end

        created_flag = true
      else
        if created_flag
          for child in contents(grids[][3]) delete!(child) end
          trim!(grids[][3])

          str = Observable(String[" ", " ", " "])
          created_flag = false
        end

        ds[] = s
      end
  end

  create_button_next_main(fig, grids[], ds, flag)

  return screen, fig, ds, rule, flag
end
