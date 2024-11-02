#import Pkg; Pkg.add("DifferentialEquations"); Pkg.add("DynamicalSystems"); Pkg.add("CairoMakie"); Pkg.add("GLMakie"); Pkg.add("IntervalTrees"); Pkg.add("ForwardDiff")
using DifferentialEquations
using DynamicalSystems
using CairoMakie
using GLMakie
using Base.Threads
using IntervalTrees
using ForwardDiff
using LinearAlgebra

function make_return_vec(v::Vector{Vector{Point2}})
    res = Vector{Vector{Vector}}()

    for i in eachindex(v)
        push!(res, [[v[i][1][2], v[i][1][1]]])
        for j in eachindex(v[i])
            if (ind = findfirst(==(res[i][j][2]), vec[2] for vec in v[i])) === nothing
                break
            else
                push!(res[i], [v[i][ind][2], v[i][ind][1]])
            end
        end
    end

    return res
end

function orbits(matrix::Vector{Vector}, t)
    s = length(t)

    vec = Vector{Vector{Point2}}()

    used = Matrix{Bool}(undef, s, s)
    fill!(used, false)

    j = 1
    while j < s
        i = j + 10
        while i < s
            if matrix[i][j]
                indices = find_orbit(matrix, i, j, s, used)
                if length(indices) > 2
                    push!(vec, indices)
                    j = next_column(indices)
                    break
                end
            end
            i += 1
        end
        j += 1
    end
    v = filter_orbits(vec)
    return sort(v)
end

function filter_orbits(vec::Vector{Vector{Point2}})
    sorted_vec = sort(vec, by = length, rev = true)
    result = Vector{Vector{Point2}}()
    treey = IntervalTree{Int, Interval{Int}}()
    treex = IntervalTree{Int, Interval{Int}}()

    for indices in sorted_vec
        if all(!hasintersection(treey, i) for i in (minimum(v[2] for v in indices)):1:(maximum(v[1] for v in indices))) &&
           all(!hasintersection(treex, v[2]) for v in indices)

            push!(result, indices)

            push!(treey, Interval(minimum(v[2] for v in indices), maximum(v[1] for v in indices)))
            push!(treex, Interval(minimum(v[2] for v in indices), maximum(v[2] for v in indices)))
        end
    end

    return result
end

function find_orbit(matrix::Vector{Vector}, i::Int64, j::Int64, s::Int64, used::Matrix{Bool})
    indices = Vector{Point2}()
    find!(matrix, i, j, indices, s, used)
    return indices
end

function find!(matrix::Vector{Vector}, i::Int64, j::Int64, indices::Vector{Point2}, s::Int64, used::Matrix{Bool})
    tmpi, tmpj = i, j
    while true
        if i < s && j < s
            if matrix[i][j] && !used[i, j]
                used[i, j] = true
                push!(indices, [i, j])
                i += 1
                j += 1
            else
                break
            end
        else
            break
        end
    end

    while true
        if 0 < tmpi && 0 < tmpj
            if matrix[tmpi][tmpj] && !used[tmpi, tmpj]
                used[tmpi, tmpj] = true
                push!(indices, [tmpi, tmpj])
                tmpi -= 1
                tmpj -= 1
            else
                break
            end
        else
            break
        end
    end

    return nothing
end

function next_column(indices::Vector{Point2})
    return maximum(v[2] for v in indices)
end

function recurrence_matrix_to_decent_matrix(R::RecurrenceMatrix, s::Int64)
    r = Vector{Vector}()
    for i in 0:s:(s * (s - 1))
        tmp = Vector{}()
        for j in (i + 1):(s + i) push!(tmp, R[j]) end
        push!(r, tmp)
    end
    return r;
end

function lorenz_rule!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
    return nothing
end

f((x, y, z)) = [p0[1] * (y - x), x * (p0[2] - z) - y, x * y - p0[3] * z]

function F(points, Y, t, i)
    return points[(i + 1) % length(points) + 1] - Y[(t[1] + i + 1) % (t[2] + 1) + 1]
end

function newton_method(points::Vector{Point3f}, Y::StateSpaceSet, times; max_iter = 1000, tolerance = 1e-1)
    u0 = points[1]
    for i in 1:max_iter

        J = ForwardDiff.jacobian(x -> f(x), [u0[1], u0[2], u0[3]])

        if det(J) == 0.0 return nothing end

        u_new = u0 - inv(J) * F(points, Y, times, i)

        if norm(u_new - points[end]) <= tolerance
            return true
        end

        u0 = u_new

    end
    return false
end

function run_simulation(choice::Int64, total_time::Float64, u0::Vector{Float64}, p0::Vector{Float64})
    if choice == 1
        lorenz = CoupledODEs(lorenz_rule!, u0, p0)
        Δt = 0.01
        Y, t = trajectory(lorenz, total_time; Ttr = 3, Δt)
    elseif choice == 2
        tspan = (0.0, total_time)
        prob = ODEProblem(lorenz_rule!, u0, tspan, p0)
        sol = solve(prob, Tsit5())
        Y = StateSpaceSet(sol.u)
        t = sol.t
    else
        error("Invalid choice. Please select 1 to use a continuous dynamic system solution or 2 to use an adaptive step solution")
    end

    R = RecurrenceMatrix(Y, 0.1)
    s = size(R)[1]
    r = recurrence_matrix_to_decent_matrix(R, s)
    v = orbits(r, t)

    return Y, t, R, v
end

function main(ARGS)
    if length(ARGS) < 2
        error("You must add two command line argument's:
               first: equal to 1 to use a continuous dynamic system solution or 2 to use an adaptive step solution
               second: total time of dynamic system evaluation")
    end

    if parse(Int64, ARGS[2]) <= 0
        error("Total time cant be less then 0")
    end

    u0 = [10, 10, 10.9]
    global p0 = [10, 28, 2.666]
    total_time = parse(Float64, ARGS[2])

    Y, t, R, v = run_simulation(parse(Int64, ARGS[1]), total_time, u0, p0)

    # points - vector, which will be used to draw found candidates for periodic orbit
    points = Vector{Point2}()
    for i in eachindex(v)
        for j in eachindex(v[i])
            push!(points, [v[i][j][2], v[i][j][1]])
        end
    end

    s = length(t)
    return_vec = make_return_vec(v)

    # Making directories and paths to them
    current_dir = pwd()
    png_path = joinpath(dirname(current_dir), "Png")
    rm(png_path, recursive=true, force=true)
    mkdir(png_path)
    mkv_path = joinpath(dirname(current_dir), "Mkv")
    rm(mkv_path, recursive=true, force=true)
    mkdir(mkv_path)

    # Visualization and animation
    CairoMakie.activate!()

    # Recurrence matrix visualization
    fig = Figure(resolution = (10000, 10000))
    ax = Axis(fig[1, 1], titlesize = 130, title = "Recurrence matrix")
    x, y = coordinates(R)
    CairoMakie.scatter!(ax, x, y; color = :black, markersize = 10)
    CairoMakie.scatter!(ax, points; color = :red, markersize = 10)
    ax.limits = ((1, s), (1, s))
    ax.aspect = 1
    save_path = joinpath(png_path, "RecurrenceMatrix.png")
    save(save_path, fig)

    GLMakie.activate!()

    # Attractor visialization and orbits animation
    fig2 = Figure(resolution=(5000, 2500))
    ax1 = Axis3(fig2[1, 1], azimuth = 0.5 * pi, titlesize = 40, title = "Original attractor")
    ax2 = Axis3(fig2[1, 2], azimuth = 0.5 * pi, titlesize = 40, title = "Reconstructed attractor")
    GLMakie.lines!(ax1, Y[:, 1], Y[:, 2], Y[:, 3]; linewidth = 1, color = :blue)

    colors = [:blue, :red, :green, :yellow, :orange, :black, :gray,
        :purple, :pink, :brown, :cyan, :magenta, :turquoise, :coral,
        :lavender, :salmon, :teal, :violet, :gold, :silver]
    k, n = 0, length(colors)
    for i in eachindex(return_vec)

        points = Vector{Vector{Point3f}}()

        for j in 1:length(return_vec[i])
            push!(points, [])
            for k in return_vec[i][j][1]:return_vec[i][j][2]
                push!(points[j], Y[k])
            end
        end

        # Chech orbit with Newton method
        j = 1
        while j <= length(points)
            nw = newton_method(points[j], Y, return_vec[i][j])
            if nw == false
                deleteat!(points, j)
                j -= 1
            elseif nw === nothing
                println("TODO: Add case handling of degenerate Jacobian")
            end
            j += 1
        end

        # No true periodic orbits left after Newton method
        if length(points) == 0 continue end

        xmin, xmax, ymin, ymax, zmin, zmax =
            minimum(p[1] for subpoints in points for p in subpoints), maximum(p[1] for subpoints in points for p in subpoints),
            minimum(p[2] for subpoints in points for p in subpoints), maximum(p[2] for subpoints in points for p in subpoints),
            minimum(p[3] for subpoints in points for p in subpoints), maximum(p[3] for subpoints in points for p in subpoints)

        # Reconstructed attractor visialization
        for j in 1:length(return_vec[i])
            k = (i + j - 1) % n + 1
            GLMakie.lines!(ax2, points[j]; linewidth = 1, color = colors[k])
        end

        # Orbit's animation
        color_ind = Observable(1)
        fig3 = Figure(resolution=(1000, 1000))
        ax3 = Axis3(fig3[1, 1])
        limits!(ax3, xmin, xmax, ymin, ymax, zmin, zmax)

        frames = 1:length(return_vec[i])
        save_path = joinpath(mkv_path, "ReconsructedAttractorOrbit$i.mp4")
        record(fig3, save_path) do io
            for frame in frames
                p_color = colors[color_ind[]]
                color_ind[] = frame % length(colors) + 1
                pointObservable = Observable(Point3f[points[frame][1]])
                for j in 2:length(points[frame])
                    pointObservable[] = push!(pointObservable[], points[frame][j])
                    GLMakie.lines!(ax3, pointObservable[]; color = p_color, markersize = 4000)
                    recordframe!(io)
                end
            end
        end
    end

    save_path = joinpath(png_path, "Attractor.png")
    save(save_path, fig2)
    return nothing
end

main(ARGS)
