include("Pages/MainPage.jl")
include("Pages/SecondPage.jl")
include("Pages/ThirdPage.jl")
include("Pages/FourthPage.jl")


function main()
  GLMakie.activate!()

  screen, fig, ds, rule, pageClosedM = initMainPage()

  display(screen, fig)

  pageClosedS = Observable(false)
  pageClosedT = Observable(false)
  total_time  = Observable(1000)
  i_orbits    = Observable(Bool[])
  u0s         = Observable(Vector{Float64}[])
  system      = Observable{Any}("")
  psos        = Observable{Any}("")

  on(pageClosedM) do _ initSecondPage(fig, ds, rule, pageClosedS, u0s, total_time, system)              end

  on(pageClosedS) do _ initThirdPage(screen, fig, u0s, total_time, system, psos, pageClosedT, i_orbits) end

  on(pageClosedT) do _ initFourthPage(screen, fig, system, psos, i_orbits, ds[], u0s[], total_time[])   end
end

main()
