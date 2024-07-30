#=
import Pkg
Pkg.add("DynamicalSystems")
Pkg.add("CairoMakie")
Pkg.add("GLMakie")
=#
using DynamicalSystems
using CairoMakie
using GLMakie

function orbits(matrix, indexmin, indexmax, tmin, tmax, s)
    vec = Vector{Vector{Point2}}()
    j = 10

    used = Matrix{Bool}(undef, s, s)
    fill!(used, false)

    while j < indexmax

        for i in j+indexmin:indexmax
            if matrix[i][j]
                indices = find_orbit(matrix, i, j, used, s)
                if length(indices) >= tmin && length(indices) <= tmax && is_original_orbit(vec, indices)
                    mysort!(indices)
                    push!(vec, indices)
                    j = next_column(indices)
                    break
                end
            end
        end

        j += 1
    end

    return vec
end

function find_orbit(matrix, i, j, used, s)
    indices = Vector{Point2}()

    push!(indices, [i, j])
    used[i, j] = true

    find!(matrix, i, j, indices, used, s)
    return indices
end

function find!(matrix, i, j, indices, used, s)
    for i1 in -1:0
        for j1 in -1:0
            if i1 == j1 == 0 || i - i1 == 0 || j - j1 == 0 || i - i1 == s || j - j1 == s
                continue
            end

            if matrix[i - i1][j - j1] && !used[i - i1, j - j1]
                push!(indices, [i - i1, j - j1])
                used[i - i1, j - j1] = true
                find!(matrix, i - i1, j - j1, indices, used, s)
            end
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

    if max < min return true end

    return false
end

function mysort!(indices)
    function quicksort(arr, low, high)
        if low < high
            pivot = partition(arr, low, high)
            quicksort(arr, low, pivot - 1)
            quicksort(arr, pivot + 1, high)
        end
    end

    function partition(arr, low, high)
        pivot = arr[high][2]
        i = low - 1
        for j in low:high-1
            if arr[j][2] <= pivot
                i += 1
                arr[i], arr[j] = arr[j], arr[i]
            end
        end
        i += 1
        arr[i], arr[high] = arr[high], arr[i]
        return i
    end

    n = length(indices)
    quicksort(indices, 1, n)
    return nothing
end

function r_matrix_to_real_matrix(r, s)
    real = Vector{Vector}()

    for i in 0:s:(s * (s - 1))
        tmp = Vector{}()
        for j in (i + 1):(s + i)
            push!(tmp, r[j])
        end
        push!(real, tmp)
    end

    return real;
end

function newton_check()

    return true
end

function lorenz_rule!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
    return nothing
end

u0 = [10, 10, 10.9]
p0 = [10, 28, 2.666]
lorenz = CoupledODEs(lorenz_rule!, u0, p0)

total_time = 200
Δt = 0.01
Ttr = 3
Y, t = trajectory(lorenz, total_time; Ttr, Δt)

R = RecurrenceMatrix(Y, 3.0)

s = size(R)[1]
r = r_matrix_to_real_matrix(R, s)

tmin = 20 / Δt # Minimum time of each orbit
tmax = 50 / Δt # Maximum time of each orbit
dt = 200       # Time close to main diagonal which we won't consider

v = orbits(r, dt, s - dt, tmin, tmax, s)

p = Vector{Point2}()
for i in eachindex(v)
    for j in eachindex(v[i])
        push!(p, [v[i][j][2], v[i][j][1]])
    end
end

CairoMakie.activate!()

current_dir = pwd()

recurrenceplot(R; ascii=true)

fig = Figure(resolution=(10000, 10000))
ax = Axis(fig[1, 1])
x, y = coordinates(R)
CairoMakie.scatter!(ax, x, y; color=:black, markersize=10)
CairoMakie.scatter!(ax, p; color=:red, markersize=10)
ax.limits = ((1, s), (1, s));
ax.aspect = 1
save_path = joinpath(current_dir, "RecurrenceMatrix.png")
save(save_path, fig)

GLMakie.activate!()

colors = [:blue, :red, :green, :yellow, :orange, :black, :white, :gray,
    :purple, :pink, :brown, :cyan, :magenta, :turquoise, :coral,
    :lavender, :salmon, :teal, :violet, :gold, :silver]
fig2 = Figure(resolution=(5000, 2500))
ax1 = Axis3(fig2[1, 1], azimuth=0.5 * pi, title="Original attractor")
GLMakie.lines!(ax1, Y[:, 1], Y[:, 2], Y[:, 3]; linewidth=1, color=:blue)
ax2 = Axis3(fig2[1, 2], azimuth=0.5 * pi, title="Reconstructed attractor")
k = 1
for i in eachindex(v)
    global k
    if k > length(colors) k = 1 end
    points = Vector{Point3f}()
    for j in 1:length(v[i])
        push!(points, [Y[v[i][j][1]][1], Y[v[i][j][1]][2], Y[v[i][j][1]][3]])
    end
    GLMakie.lines!(ax2, points; linewidth=1, color=colors[k])
    k += 1
end
save_path = joinpath(current_dir, "Attractor.png")
save(save_path, fig2)

mkv_path = joinpath(current_dir, "Mkv")
rm(mkv_path, recursive=true, force=true)
mkdir(mkv_path)

for i in eachindex(v)
    points = Vector{Point3f}()
    for j in eachindex(v[i])
        push!(points, [Y[v[i][j][1]][1], Y[v[i][j][1]][2], Y[v[i][j][1]][3]])
    end

    xmin, xmax, ymin, ymax, zmin, zmax =
        minimum(p[:][1] for p in points), maximum(p[:][1] for p in points),
        minimum(p[:][2] for p in points), maximum(p[:][2] for p in points),
        minimum(p[:][3] for p in points), maximum(p[:][3] for p in points)

    point = Observable(Point3f[points[1]])
    figi = Figure(resolution=(1000, 1000))
    axi = Axis3(figi[1, 1], title="$i orbit")
    limits!(axi, xmin, xmax, ymin, ymax, zmin, zmax)
    GLMakie.lines!(axi, point, color=:blue, markersize = 4000)
    
    frames = 2:length(points)
    save_path = joinpath(mkv_path, "ReconsructedAttractorOrbit$i.mp4")
    record(figi, save_path, frames;) do frame
        point[] = push!(point[], points[frame])
    end
end
