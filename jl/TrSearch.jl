#import Pkg; Pkg.add("DynamicalSystems"); Pkg.add("CairoMakie"); Pkg.add("GLMakie"); Pkg.add("IntervalTrees")
using DynamicalSystems
using CairoMakie
using GLMakie
using Base.Threads
using IntervalTrees

function orbits(matrix, indexmin, tmin, tmax, s)
    indexmax = s - indexmin

    vec = Vector{Vector{Point2}}()

    used = Matrix{Bool}(undef, s, s)
    fill!(used, false)

    local_vecs = [Vector{Any}() for _ in 1:nthreads()]

    Threads.@threads for j in indexmin:indexmax
        for i in j + indexmin:indexmax
            if matrix[i][j]
                indices = find_orbit(matrix, i, j, s, used)
                if length(indices) >= tmin && length(indices) <= tmax
                    push!(local_vecs[threadid()], indices)
                    j = next_column(indices)
                    break
                end
            end
        end
    end

    for local_vec in local_vecs
        append!(vec, local_vec)
    end

    return filter_orbits(vec)
end

function filter_orbits(vec)
    sorted_vec = sort(vec, by = length, rev = true)
    result = Vector{Vector{Point2}}()
    treex = IntervalTree{Int, Interval{Int}}()
    treey = IntervalTree{Int, Interval{Int}}()

    for indices in sorted_vec
        if all(!hasintersection(treex, v[2]) for v in indices) &&
            all(!hasintersection(treey, v[1]) for v in indices)

            push!(result, indices)

            push!(treex, Interval(minimum(v[2] for v in indices), maximum(v[2] for v in indices)))
            push!(treey, Interval(minimum(v[1] for v in indices), maximum(v[1] for v in indices)))
        end
    end
        
    return result
end

function find_orbit(matrix, i, j, s, used)
    indices = Vector{Point2}()
    find!(matrix, i, j, indices, s, used)
    return indices
end

function find!(matrix, i, j, indices, s, used)
    if 0 < i < s && 0 < j < s
        if matrix[i][j] && !used[i, j]
            used[i, j] = true
            push!(indices, [i, j])
            find!(matrix, i + 1, j + 1, indices, s, used)
            find!(matrix, i - 1, j - 1, indices, s, used)
        end
    end
    return nothing
end

function next_column(indices)
    return maximum(v[2] for v in indices)
end

function is_original_orbit(vec, indices)
    n = length(vec)
    if n == 0 return true end

    max = maximum(v[:][2] for v in vec[n])
    min = minimum(v[2] for v in indices)

    if max < min return true
    else return false
    end
end

function r_matrix_to_real_matrix(r, s)
    real = Vector{Vector}()
    for i in 0:s:(s * (s - 1))
        tmp = Vector{}()
        for j in (i + 1):(s + i) push!(tmp, r[j]) end
        push!(real, tmp)
    end
    return real;
end

function lorenz_rule!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
    return nothing
end

function main()
    u0 = [10, 10, 10.9]
    p0 = [10, 28, 2.666]
    lorenz = CoupledODEs(lorenz_rule!, u0, p0)

    Δt = 0.01
    total_time = 200
    Y, t = trajectory(lorenz, total_time; Ttr = 3, Δt)

    R = RecurrenceMatrix(Y, 3.0)

    s = size(R)[1]
    r = r_matrix_to_real_matrix(R, s)

    tmin = 2 / Δt  # Minimum time of each orbit
    tmax = 40 / Δt # Maximum time of each orbit
    dt = 200       # Time close to main diagonal which we won't consider

    # Main alorithm
    v = orbits(r, dt, tmin, tmax, s)

    points = Vector{Point2}()
    for i in eachindex(v)
        for j in eachindex(v[i])
            push!(points, [v[i][j][2], v[i][j][1]])
        end
    end

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

    fig2 = Figure(resolution=(7500, 2500))
    ax1 = Axis3(fig2[1, 1], azimuth = 0.5 * pi, titlesize = 40, title = "Original attractor")
    ax2 = Axis3(fig2[1, 2], azimuth = 0.5 * pi, titlesize = 40, title = "Reconstructed attractor 1")
    ax3 = Axis3(fig2[1, 3], azimuth = 0.5 * pi, titlesize = 40, title = "Reconstructed attractor 2")
    GLMakie.lines!(ax1, Y[:, 1], Y[:, 2], Y[:, 3]; linewidth = 1, color = :blue)

    colors = [:blue, :red, :green, :yellow, :orange, :black, :gray,
        :purple, :pink, :brown, :cyan, :magenta, :turquoise, :coral,
        :lavender, :salmon, :teal, :violet, :gold, :silver]
    k, n = 0, length(colors)
    for i in eachindex(v)
        k += 1
        if k > n k = 1 end

        pointsi, pointsj, timei, timej = 
            Vector{Point3f}(), Vector{Point3f}(), Vector{}(), Vector{}()

        for j in 1:length(v[i])
            push!(timei, v[i][j][1])
            push!(timej, v[i][j][2])
        end
        sort!(timei)
        sort!(timej)

        for j in 1:length(v[i])
            push!(pointsi, [Y[timei[j]][1], Y[timei[j]][2], Y[timei[j]][3]])
            push!(pointsj, [Y[timej[j]][1], Y[timej[j]][2], Y[timej[j]][3]])
        end

        GLMakie.lines!(ax2, pointsi; linewidth = 1, color = colors[k])
        GLMakie.lines!(ax3, pointsj; linewidth = 1, color = colors[k])

        xmin, xmax, ymin, ymax, zmin, zmax =
            minimum(p[:][1] for p in pointsi), maximum(p[:][1] for p in pointsi),
            minimum(p[:][2] for p in pointsi), maximum(p[:][2] for p in pointsi),
            minimum(p[:][3] for p in pointsi), maximum(p[:][3] for p in pointsi)

        xmin1, xmax1, ymin1, ymax1, zmin1, zmax1 =
            minimum(p[:][1] for p in pointsj), maximum(p[:][1] for p in pointsj),
            minimum(p[:][2] for p in pointsj), maximum(p[:][2] for p in pointsj),
            minimum(p[:][3] for p in pointsj), maximum(p[:][3] for p in pointsj)
        
        point1 = Observable(Point3f[pointsi[1]])
        point2 = Observable(Point3f[pointsj[1]])
        figi = Figure(resolution=(1000, 1000))
        axi = Axis3(figi[1, 1])
        axj = Axis3(figi[1, 2])
        limits!(axi, xmin, xmax, ymin, ymax, zmin, zmax)
        limits!(axj, xmin1, xmax1, ymin1, ymax1, zmin1, zmax1)
        
        GLMakie.lines!(axi, point1, color=:blue, markersize = 4000)
        GLMakie.lines!(axj, point2, color=:red, markersize = 4000)

        frames = 2:length(pointsi)
        save_path = joinpath(mkv_path, "ReconsructedAttractorOrbit$i.mp4")
        record(figi, save_path, frames) do frame
            point1[] = push!(point1[], pointsi[frame])
            point2[] = push!(point2[], pointsj[frame])
        end
    end

    save_path = joinpath(png_path, "Attractor.png")
    save(save_path, fig2)
end

main()