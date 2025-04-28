using DynamicalSystems, OrdinaryDiffEq
using GLMakie, WGLMakie, Makie.Colors
using GridLayoutBase

include("util.jl")


function convert_to_svector(points::Vector{Point3f})
    return [SVector(Float64(p[1]), Float64(p[2]), Float64(p[3])) for p in points]
end

function scan_poincaresos(
    As::Vector{Tuple{S, T}} where {S<:AbstractStateSpaceSet, T<:StepRangeLen},
    args...; kwargs...)
    return scan_poincaresos([As[i][1] for i in eachindex(As)], args...; kwargs...)
end

function scan_poincaresos(
    As::Vector{<:AbstractStateSpaceSet}, screen, current_size, best_psos, i_orbits, total_time, pageClosedT;
    linekw = (), scatterkw = (), direction = -1, j = 2,
    colors = [:red, :orange, :green, :yellow, :blue, :pink,
              :brown, :cyan, :magenta, :purple, :coral, :gold,
              :lavender, :salmon, :teal, :gray, :violet, :silver])
    
    size_ = size(As, 1)
    @assert size(colors, 1) > size_

    was_changed      = false
    psos_active      = false
    orbit_number     = nothing
    buttonPsos       = nothing
    buttonChangePsos = nothing

    mi, ma, minmax = total_minmaxima(As)

    otheridxs = [setdiff(1:3, j)...]

    figure = Figure(size = current_size)
    display(screen, figure)

    ax  = Axis3(figure[1, 1])
    axp = Axis(figure[1, 2])

    GLMakie.limits!(ax,  mi[1], ma[1], mi[2], ma[2], mi[3], ma[3])
    GLMakie.limits!(axp, mi[otheridxs[1]], ma[otheridxs[1]],
                         mi[otheridxs[2]], ma[otheridxs[2]])

    sg = SliderGrid(figure[2, :],
        (label = "$(('x':'z')[j]) =", range = range(mi[j], ma[j]; length = total_time),
        startvalue = (ma[j] + mi[j]) / 2)
    )
    slider = sg.sliders[1]

    line = []
    @inbounds for i in eachindex(As)
        push!(line, lines!(ax, As[i].data; color = colors[i], transparency = false, linekw...))
    end

    scatter_plots_ax  = Observable(Point3f[])
    scatter_plots_axp = Observable(Point2f[])
    point_color       = Observable(Symbol[])

    axp_scatter = Makie.scatter!(axp, scatter_plots_axp; color = point_color, scatterkw...)
    ax_scatter  = Makie.scatter!(ax, scatter_plots_ax; color = point_color, markersize = 5, scatterkw...)

    orbits_range = eachindex(As)

    # Main animation
    on(slider.value) do y_val
        all_psos3d = Point3f[]
        all_psos2d = Point2f[]
        all_colors = Symbol[]
        @inbounds for i in orbits_range
            psos = nothing

            try
                psos = poincaresos(As[i], (j, y_val); direction, warning = false)
            catch err
                continue
            end

            isnothing(psos) && continue

            filtered = filter(!isnothing, psos)

            if !isempty(filtered)
                psos2d = [p[otheridxs] for p in filtered]
                psos3d = [p for p in filtered]

                append!(all_psos2d, psos2d)
                append!(all_psos3d, psos3d)

                @inbounds for _ in filtered
                    push!(all_colors, colors[i])
                end
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

    # Buttons
    if (size_ != 1)
        buttongrid = GridLayout(figure[3, :])

        # Button to reset everything
        buttonResetAll = Button(buttongrid[1, 1], label = "Reset all")
        on(buttonResetAll.clicks) do _
            GLMakie.limits!(ax,  mi[1], ma[1], mi[2], ma[2], mi[3], ma[3])
            GLMakie.limits!(axp, mi[otheridxs[1]], ma[otheridxs[1]],
                                 mi[otheridxs[2]], ma[otheridxs[2]])

            for k in eachindex(line)
                delete!(screen, figure.scene, line[k])
            end
            empty!(line)

            psos_active == true && (delete!(buttonPsos); delete!(buttonChangePsos);
                                    deletecol!(buttongrid, size_ + 3); deletecol!(buttongrid, size_ + 4);
                                    psos_active = false)

            delete!(screen, figure.scene, ax_scatter)
            delete!(screen, figure.scene, axp_scatter)

            scatter_plots_ax  = Observable(Point3f[])
            scatter_plots_axp = Observable(Point2f[])
            point_color       = Observable(Symbol[])

            @inbounds for i in eachindex(As)
                push!(line, lines!(ax, As[i].data; color = colors[i], transparency = false, linekw...))
            end

            ss = copy(mi)
            ws = copy(ma) .- copy(mi)
            ws[j] = 0

            was_changed == true && (was_changed = false)

            slider.range = range(mi[j], ma[j]; length = total_time)
            set_close_to!(slider, (mi[j] + ma[j]) / 2)
            
            was_changed = true

            orbits_range = eachindex(As)
            orbit_number = nothing
            ss = copy(mi)
            ws = copy(ma) .- copy(mi)
            ws[j] = 0
        end

        # Buttons to change orbits
        vGrids = [vbox!(Button(figure, label = "Change to $i"), Checkbox(figure, checked = true)) for i in eachindex(As)]
        for i in eachindex(As)
            buttongrid[1, 1 + i] = vGrids[i]
            on(vGrids[i].content[1].content.clicks) do _
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

                scatter_plots_ax  = Observable(Point3f[])
                scatter_plots_axp = Observable(Point2f[])
                point_color       = Observable(Symbol[])

                push!(line, lines!(ax, As[i].data; color = :turquoise, transparency = false, linekw...))

                ss = copy(minmax[i][1])
                ws = copy(minmax[i][2]) .- copy(minmax[i][1])
                ws[j] = 0

                was_changed == true && (was_changed = false)

                slider.range = range(minmax[i][1][j], minmax[i][2][j]; length = total_time)
                set_close_to!(slider, (minmax[i][1][j] + minmax[i][2][j]) / 2)
                was_changed = true

                orbits_range = range(i, i)
                orbit_number = i

                if !psos_active
                    GridLayoutBase.insertcols!(buttongrid, size_ + 3, 2)
                    # Button to change optimal psos
                    buttonChangePsos = Button(buttongrid[1, size_ + 2], label = "Change optimal psos")
                    on(buttonChangePsos.clicks) do _
                        length(scatter_plots_ax[]) != 0 &&
                            (best_psos[][orbit_number][j] = 
                                (s = StateSpaceSet(convert_to_svector(scatter_plots_ax[])),
                                 v = scatter_plots_ax[][1][j]))
                    end

                    # Button to show optimal best_psos
                    buttonPsos = Button(buttongrid[1, size_ + 3], label = "Show optimal psos")
                    on(buttonPsos.clicks) do _ set_close_to!(slider, best_psos[][orbit_number][j].v) end
                    psos_active = true
                end
            end

            # Checkboxes
            on(vGrids[i].content[2].content.checked) do b
              i_orbits[][i] = b
            end 
        end

        # Button to change direction
        buttonChangeJ = Button(buttongrid[1, size_ + 4], label = "Change poincaresos direction")
        on(buttonChangeJ.clicks) do _
            j = j + 1
            j == 4 && (j = 1)
            otheridxs = [setdiff(1:3, j)...]
            if isnothing(orbit_number)
                GLMakie.limits!(axp, mi[otheridxs[1]], ma[otheridxs[1]],
                                     mi[otheridxs[2]], ma[otheridxs[2]])

                ss = copy(mi)
                ws = copy(ma) .- copy(mi)
                ws[j] = 0 # Should be changes to (min + max) / 2
                delete!(screen, figure.scene, ax_scatter)
                delete!(screen, figure.scene, axp_scatter)

                scatter_plots_ax  = Observable(Point3f[])
                scatter_plots_axp = Observable(Point2f[])
                point_color       = Observable(Symbol[])

                was_changed == true && (was_changed = false)

                slider.range = range(mi[j], ma[j]; length = total_time)
                set_close_to!(slider, (mi[j] + ma[j]) / 2)
                was_changed = true
            else
                ss = copy(minmax[orbit_number][1])
                ws = copy(minmax[orbit_number][2]) .- copy(minmax[orbit_number][1])
                ws[j] = 0 # Should be changes to (min + max) / 2

                GLMakie.limits!(axp, minmax[orbit_number][1][otheridxs[1]], minmax[orbit_number][2][otheridxs[1]],
                                     minmax[orbit_number][1][otheridxs[2]], minmax[orbit_number][2][otheridxs[2]])

                delete!(screen, figure.scene, ax_scatter)
                delete!(screen, figure.scene, axp_scatter)

                scatter_plots_ax  = Observable(Point3f[])
                scatter_plots_axp = Observable(Point2f[])
                point_color       = Observable(Symbol[])

                was_changed == true && (was_changed = false)

                slider.range = range(minmax[orbit_number][1][j], minmax[orbit_number][2][j]; length = total_time)
                set_close_to!(slider, (minmax[orbit_number][1][j] + minmax[orbit_number][2][j]) / 2)
                was_changed = true
            end
        end
    end
    
    # Button to go to the next page
    buttonNext = Button(buttongrid[1, size_ + 5], label = "Next")
    on(buttonNext.clicks) do _
      delete_all!(figure)
      pageClosedT[] = true
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

    return figure, best_psos
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


function find_best_poincaresos(As::Vector{Tuple{S, T}} where {S<:AbstractStateSpaceSet, T<:StepRangeLen}, D::Int64, total_time)
    return find_best_poincaresos([As[i][1] for i in eachindex(As)], D, total_time)
end

function find_best_poincaresos(As::Vector{<:AbstractStateSpaceSet}, D::Int64, total_time)
    psos = []
    for i in eachindex(As)
        psos_j = Vector{Any}(undef, D)
 
        min, max = minmaxima(As[i])

        for j in 1:D
            psos_j[j] = (s = AbstractStateSpaceSet{D, Float64}[], v = Float64)
            
            # TODO: what step should we use?
            # benc: t_t = 1000, lorenz, 3 u0s, dt = 0.01
            # 119s - all dots
            # 42s  - all dots + threads
            # 5s   - fixed step, Δ = 0.01
            # 2.7s - fixed step, Δ = 0.01 + threads
            # 1.2s - "adaptive" step
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
        push!(psos, psos_j)
    end

    return psos
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
