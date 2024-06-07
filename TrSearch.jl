#import Pkg; Pkg.add("DynamicalSystems"); Pkg.add("CairoMakie"); Pkg.add("GLMakie")
using DynamicalSystems
using CairoMakie
using GLMakie

function orbits(matrix, indexmin, indexmax, tmin, tmax, s)
    vec = Vector{Vector{Vector}}()
    j = 10

    used = Matrix{Bool}(undef, s, s)
    fill!(used, false)

    while true
        if j >= indexmax
            break
        end

        for i in j + indexmin:indexmax
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
    indices = Vector{Vector}()

    push!(indices, [i, j])
    used[i, j] = true

    find!(matrix, i, j, indices, used, s)
    return indices
end

function find!(matrix, i, j, indices, used, s)
    for i1 in -1:1
        for j1 in -1:1
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
    if n == 0
        return true
    end

    max = maximum(v[:][2] for v in vec[n])

    min = minimum(v[2] for v in indices)

    if max < min
        return true
    end
    return false
end

function mysort!(indices)
    n = length(indices)
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
        for j in low:high - 1
            if arr[j][2] <= pivot
                i += 1
                arr[i][1], arr[j][1] = arr[j][1], arr[i][1]
                arr[i][2], arr[j][2] = arr[j][2], arr[i][2]
            end
        end
        i += 1
        arr[i][1], arr[high][1] = arr[high][1], arr[i][1]
        arr[i][2], arr[high][2] = arr[high][2], arr[i][2]
        return i
    end

    quicksort(indices, 1, n)
    return nothing
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

r = Vector{Vector}()
s = size(R)[1]
for i in 0:s:(s*(s - 1))
   tmp = Vector{}()
   for j in (i + 1):(s + i)
       push!(tmp, R[j])
   end
   push!(r, tmp)
end

tmin = 30 / Δt # Minimum time of each orbit
tmax = 50 / Δt # Maximum time of each orbit
dt = 200       # Time close to main diagonal which we wont consider

v = orbits(r, dt, s - dt, tmin, tmax, s)

x1 = Vector{Int64}()
y1 = Vector{Int64}()
for i in eachindex(v)
   for j in eachindex(v[i])
       push!(x1, v[i][j][2])
       push!(y1, v[i][j][1])
   end
end

CairoMakie.activate!()

current_dir = pwd()

recurrenceplot(R; ascii = true)

fig = Figure(resolution = (10000, 10000))
ax = Axis(fig[1,1])
x, y = coordinates(R)
CairoMakie.scatter!(ax, x, y; color = :black, markersize = 10)
CairoMakie.scatter!(ax, x1, y1; color = :red, markersize = 10)
ax.limits = ((1, s), (1, s));
ax.aspect = 1
save_path = joinpath(current_dir, "RecurrenceMatrix.png")
save(save_path, fig)

GLMakie.activate!()

x = Float64[]
y = Float64[]
z = Float64[]
for i in eachindex(v)
    for j in 1:length(v[i])
        push!(x, Y[v[i][j][1]][1])
        push!(y, Y[v[i][j][1]][2])
        push!(z, Y[v[i][j][1]][3])
    end
end

fig2 = Figure(resolution = (1000, 500))
ax1 = Axis3(fig2[1, 1], azimuth = 0.5 * pi, title = "Original attractor")
GLMakie.lines!(ax1, Y[:, 1], Y[:, 2], Y[:, 3]; linewidth = 1, color = :blue)
ax2 = Axis3(fig2[1, 2], azimuth = 0.5 * pi, title = "Reconstructed attractor")
GLMakie.lines!(ax2, x, y, z; linewidth = 1, color = :red)
save_path = joinpath(current_dir, "Attractor.png")
save(save_path, fig2)

png_path = joinpath(current_dir, "Png")
rm(png_path, recursive = true, force = true)
mkdir(png_path)

for i in eachindex(v)
    x = Float64[]
    y = Float64[]
    z = Float64[]
    for j in eachindex(v[i])
        push!(x, Y[v[i][j][1]][1])
        push!(y, Y[v[i][j][1]][2])
        push!(z, Y[v[i][j][1]][3])
    end

    figi = Figure(resolution = (1000, 1000))
    axi = Axis3(figi[1, 1], title = "$i orbit")
    GLMakie.lines!(axi, x, y, z; linewidth = 5, color = :blue)
    save_path = joinpath(current_dir, "Png\\ReconsructedAttractorOrbit$i.png")
    save(save_path, figi)
end