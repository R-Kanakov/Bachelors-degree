using DynamicalSystems, OrdinaryDiffEq
using GLMakie, WGLMakie, Makie.Colors
using GridLayoutBase


function convert_to_svector(points::Vector{Point3f})
    return [SVector(Float64(p[1]), Float64(p[2]), Float64(p[3])) for p in points]
end

function scan_poincaresos(
    As::Vector{Tuple{S, T}} where {S<:AbstractStateSpaceSet, T<:StepRangeLen},
    args...; kwargs...
)
    return scan_poincaresos([As[i][1] for i in eachindex(As)], args...; kwargs...)
end
# Correct dir when final psos shown and change it
function scan_poincaresos(
    As::Vector{<:AbstractStateSpaceSet},
    screen, current_size, best_psos, final_psos, i_orbits, total_time, pageClosedT;
    linekw = (), scatterkw = (), j = 2,
    colors = [:red, :orange, :green, :yellow, :pink, :brown,
              :cyan, :magenta, :purple, :coral, :gold, :teal,
              :salmon, :lavender, :gray, :violet, :silver]
)
    size_ = size(As, 1)
    @assert size(colors, 1) > size_

    each_i             = eachindex(As)
    was_changed        = false
    orbit_number       = nothing
    button_change_psos = nothing

    mi, ma, minmax = total_minmaxima(As)

    otheridxs = [setdiff(1:3, j)...]

    figure = Figure(size = current_size)

    display(screen, figure)

    ax  = Axis3(figure[1, 1])
    axp = Axis(figure[1, 2])

    GLMakie.limits!(ax,  mi[1], ma[1],
                         mi[2], ma[2],
                         mi[3], ma[3])
    GLMakie.limits!(axp, mi[otheridxs[1]], ma[otheridxs[1]],
                         mi[otheridxs[2]], ma[otheridxs[2]])

    sg = SliderGrid(figure[2, :],
        (label     = "$(('x':'z')[j]) =",
        range      = range(mi[j], ma[j], length = total_time),
        startvalue = (ma[j] + mi[j]) / 2)
    )
    slider = sg.sliders[1]

    if size_ == 1
        lines!(ax, As[1].data; color = RGBA{Float32}(0.0, 0.44705883, 0.69803923, 1.0), transparency = false, linekw...)
    else
        line = []
        @inbounds for i in each_i
            push!(line, lines!(ax, As[i].data; color = colors[i], transparency = false, linekw...))
        end
    end

    direction         = Observable(-1)
    scatter_plots_ax  = Observable(Point3f[])
    scatter_plots_axp = Observable(Point2f[])
    point_color       = Observable(Symbol[])
    curr_plane        = Observable{Any}()

    final_psos[] = (psos = [best_psos[k][j].s for k in eachindex(best_psos)], v = best_psos[1][j].v, j = j, d = direction[])

    axp_scatter = Makie.scatter!(axp, scatter_plots_axp; color = point_color, scatterkw...)
    ax_scatter  = Makie.scatter!(ax,  scatter_plots_ax;  color = point_color, markersize = 5, scatterkw...)

    orbits_range = each_i

    # Main animation
    on(slider.value) do j_val
        all_psos3d = Point3f[]
        all_psos2d = Point2f[]
        all_colors = Symbol[]

        curr_plane[] = (j, j_val)

        @inbounds for i in orbits_range
            psos = nothing

            try
                psos = poincaresos(As[i], curr_plane[]; direction = direction[], warning = false)
            catch err
                continue
            end

            psos2d = [p[otheridxs] for p in psos]
            append!(all_psos2d, psos2d)

            psos3d = [p for p in psos]
            append!(all_psos3d, psos3d)

            @inbounds for _ in psos
                push!(all_colors, colors[i])
            end
        end

        scatter_plots_axp[] = all_psos2d
        scatter_plots_ax[]  = all_psos3d
        point_color[]       = all_colors

        if was_changed
            axp_scatter = Makie.scatter!(axp, scatter_plots_axp; color = point_color, scatterkw...)
            ax_scatter  = Makie.scatter!(ax,  scatter_plots_ax;  color = point_color, scatterkw..., markersize = 5)
            was_changed = false
        end
    end

    # Surface animation
    ss = copy(mi)
    ws = copy(ma) .- copy(mi)
    ws[j] = 0

    p = map(slider.value) do val
        ss[j] = val
        o = Makie.Point3f(ss...)
        w = Makie.Point3f(ws...)
        Makie.Rect3f(o, w)
    end

    mesh!(ax, p; color = RGBAf(0.2, 0.2, 0.25, 0.5), transparency = true)

    # Buttons
    buttons_grid = GridLayout(figure[3, :])

    if size_ != 1
        # Button to reset everything
        button_reset = Button(buttons_grid[1, 1], label = "Reset all")
        on(button_reset.clicks) do _
            GLMakie.limits!(ax,  mi[1], ma[1],
                                 mi[2], ma[2],
                                 mi[3], ma[3])
            GLMakie.limits!(axp, mi[otheridxs[1]], ma[otheridxs[1]],
                                 mi[otheridxs[2]], ma[otheridxs[2]])

            for k in eachindex(line)
                delete!(screen, figure.scene, line[k])
            end
            empty!(line)

            delete!(screen, figure.scene, ax_scatter)
            delete!(screen, figure.scene, axp_scatter)

            @inbounds for i in each_i
                push!(line, lines!(ax, As[i].data; color = colors[i], transparency = false, linekw...))
            end

            was_changed == true && (was_changed = false)

            slider.range = range(mi[j], ma[j]; length = total_time)
            set_close_to!(slider, (mi[j] + ma[j]) / 2)

            was_changed = true

            scatter_plots_ax  = Observable(Point3f[])
            scatter_plots_axp = Observable(Point2f[])
            point_color       = Observable(Symbol[])

            orbits_range = each_i
            orbit_number = nothing

            ss = copy(mi)
            ws = copy(ma) .- copy(mi)
            ws[j] = 0
        end

        # Buttons to change orbits
        vertical_grids = [vbox!(Button(figure, label = "Change to $i"),
                                Checkbox(figure, checked = true))
                          for i in each_i]

        @inbounds for i in each_i
            buttons_grid[1, 1 + i] = vertical_grids[i]
            on(vertical_grids[i].content[1].content.clicks) do _
                GLMakie.limits!(ax,  minmax[i][1][1], minmax[i][2][1], minmax[i][1][2],
                                     minmax[i][2][2], minmax[i][1][3], minmax[i][2][3])
                GLMakie.limits!(axp, minmax[i][1][otheridxs[1]], minmax[i][2][otheridxs[1]],
                                     minmax[i][1][otheridxs[2]], minmax[i][2][otheridxs[2]])

                for k in eachindex(line)
                    delete!(screen, figure.scene, line[k])
                end
                empty!(line)

                delete!(screen, figure.scene, ax_scatter)
                delete!(screen, figure.scene, axp_scatter)

                push!(line, lines!(ax, As[i].data; color = RGBA{Float32}(0.0, 0.44705883, 0.69803923, 1.0), transparency = false, linekw...))

                ss = copy(minmax[i][1])
                ws = copy(minmax[i][2]) .- copy(minmax[i][1])
                ws[j] = 0

                was_changed == true && (was_changed = false)

                slider.range = range(minmax[i][1][j], minmax[i][2][j]; length = total_time)
                set_close_to!(slider, (minmax[i][1][j] + minmax[i][2][j]) / 2)

                was_changed = true

                scatter_plots_ax  = Observable(Point3f[])
                scatter_plots_axp = Observable(Point2f[])
                point_color       = Observable(Symbol[])

                orbits_range = range(i, i)
                orbit_number = i
            end

            # Checkboxes
            on(vertical_grids[i].content[2].content.checked) do b
              i_orbits[][i] = b
            end 
        end

        # Button to change final psos
        button_change_psos = Button(buttons_grid[1, size_ + 2], label = "Change final psos")
        on(button_change_psos.clicks) do _
            if length(scatter_plots_ax[]) != 0
                psos_tmp = []
                for i in eachindex(As)
                    try
                        push!(psos_tmp, poincaresos(As[i], curr_plane[]; direction = direction[], warning = false))
                    catch err end
                end
                final_psos[] = (psos = psos_tmp, v = slider.value[], j = j, d = direction[])
            end
        end

        # Button to show final psos
        button_change_psos = Button(buttons_grid[1, size_ + 3], label = "Show final psos")
        on(button_change_psos.clicks) do _
            j != final_psos[].j && (button_change_hyperplane.clicks[] = button_change_hyperplane.clicks[] + 1)
            j != final_psos[].j && (button_change_hyperplane.clicks[] = button_change_hyperplane.clicks[] + 1)
            set_close_to!(slider, final_psos[].v)
        end
    else
        # Button to change final psos
        button_change_psos = Button(buttons_grid[1, 1], label = "Change final psos")
        on(button_change_psos.clicks) do _
            length(scatter_plots_ax[]) != 0 && 
                (final_psos[] = (psos = scatter_plots_ax[], v = slider.value[], j = j, d = direction[]))
        end

        # Button to show final psos
        button_change_psos = Button(buttons_grid[1, 2], label = "Show final psos")
        on(button_change_psos.clicks) do _
            j != final_psos[].j && (button_change_hyperplane.clicks[] = button_change_hyperplane.clicks[] + 1)
            j != final_psos[].j && (button_change_hyperplane.clicks[] = button_change_hyperplane.clicks[] + 1)
            set_close_to!(slider, final_psos[].v)
        end
    end

    # Button to change hyperplane direction
    button_change_hyperplane = Button(buttons_grid[1, size_ + 4], label = "Change hyperplane direction")
    on(button_change_hyperplane.clicks) do _
        j = (j + 1 == 4) ? 1 : j + 1
        otheridxs = [setdiff(1:3, j)...]
        if isnothing(orbit_number)
            GLMakie.limits!(axp, mi[otheridxs[1]],
                                 ma[otheridxs[1]],
                                 mi[otheridxs[2]],
                                 ma[otheridxs[2]])

            ss = copy(mi)
            ws = copy(ma) .- copy(mi)
            ws[j] = 0

            delete!(screen, figure.scene, ax_scatter)
            delete!(screen, figure.scene, axp_scatter)

            was_changed == true && (was_changed = false)

            slider.range = range(mi[j], ma[j]; length = total_time)
            set_close_to!(slider, (mi[j] + ma[j]) / 2)

            scatter_plots_ax  = Observable(Point3f[])
            scatter_plots_axp = Observable(Point2f[])
            point_color       = Observable(Symbol[])

            was_changed = true
        else
            ss = copy(minmax[orbit_number][1])
            ws = copy(minmax[orbit_number][2]) .- copy(minmax[orbit_number][1])
            ws[j] = 0

            GLMakie.limits!(axp, minmax[orbit_number][1][otheridxs[1]],
                                 minmax[orbit_number][2][otheridxs[1]],
                                 minmax[orbit_number][1][otheridxs[2]],
                                 minmax[orbit_number][2][otheridxs[2]])

            delete!(screen, figure.scene, ax_scatter)
            delete!(screen, figure.scene, axp_scatter)

            was_changed == true && (was_changed = false)

            slider.range = range(minmax[orbit_number][1][j], minmax[orbit_number][2][j]; length = total_time)
            set_close_to!(slider, (minmax[orbit_number][1][j] + minmax[orbit_number][2][j]) / 2)

            scatter_plots_ax  = Observable(Point3f[])
            scatter_plots_axp = Observable(Point2f[])
            point_color       = Observable(Symbol[])

            was_changed = true
        end
    end

    # Button to change psos direction
    button_change_direction = Button(buttons_grid[1, size_ + 5], label = "Change poincaresos direction")
    on(button_change_direction.clicks) do _
        if direction[] == 1
            direction[] = -1
        elseif direction[] == -1
            direction[] = 1
        end

        if length(scatter_plots_ax[]) != 0
            tmp = slider.value[]
            slider.value[] = 0
            set_close_to!(slider, tmp)
        end
    end

    # Button to go to the next page
    button_next = Button(buttons_grid[1, size_ + 6], label = "Next")
    on(button_next.clicks) do _
      delete_all!(figure)
      pageClosedT[] = true
    end

    return best_psos
end


function Base.filter(f, t::StateSpaceSet)
    filtered_data = [point for point in t if f(point)]
    return StateSpaceSet(filtered_data)
end

function total_minmaxima(As::Vector{<:AbstractStateSpaceSet})
    minmax = Vector{}()
    mi, ma = DynamicalSystems.minmaxima(As[1])

    mi = Vector(mi); ma = Vector(ma)

    push!(minmax, [mi, ma])

    @inbounds for j in 2:length(As)
        mi2, ma2 = DynamicalSystems.minmaxima(As[j])
        mi2 = Vector(mi2)
        ma2 = Vector(ma2)

        mi = min.(mi, mi2)
        ma = max.(ma, ma2)
    
        push!(minmax, [mi2, ma2])
    end

    return mi, ma, minmax
end


function find_best_poincaresos(As::Vector{Tuple{S, T}} where {S<:AbstractStateSpaceSet, T<:StepRangeLen}, D::Int, total_time)
    return find_best_poincaresos([As[i][1] for i in eachindex(As)], D, total_time)
end

function find_best_poincaresos(As::Vector{<:AbstractStateSpaceSet}, D::Int, total_time)
    best_psos = []
    for i in eachindex(As)
        psos_j = Vector{Any}(undef, D)

        min, max = minmaxima(As[i])

        for j in 1:D
            psos_j[j] = (s = AbstractStateSpaceSet{D, Float64}[], v = Float64)
            for j_val in range(min[j], max[j], length = total_time)
                section = nothing

                try
                    section = poincaresos(As[i], (j, j_val); direction = -1, warning = false)
                catch err
                    continue
                end

                if isempty(psos_j[j].s) || length(section) > length(psos_j[j].s)
                    psos_j[j] = (s = section, v = j_val)
                end
            end
        end
        push!(best_psos, psos_j)
    end

    return best_psos
end


function filter_orbits(trs::Vector{Tuple{S, T}} where {S<:AbstractStateSpaceSet, T<:StepRangeLen})
    return filter_orbits([trs[i][1] for i in eachindex(trs)])
end

function filter_orbits(trs::Vector{<:AbstractStateSpaceSet})
    filtered_trs = Vector{StateSpaceSet{}}()
    for tr in trs
        invalid_index = findfirst(tr) do state
            !(all(-10e6 .< state .< 10e6) && all(!isinf, state) && all(!isnan, state))
        end

        if invalid_index !== nothing
            valid_points = tr[1:invalid_index-1]
        else
            valid_points = tr
        end

        push!(filtered_trs, valid_points)
    end
    return filtered_trs
end
