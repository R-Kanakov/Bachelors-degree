#import Pkg; Pkg.add("DynamicalSystems"); Pkg.add("CairoMakie"); Pkg.add("GLMakie")
using DynamicalSystems
using CairoMakie
using GLMakie

function orbits(matrix, indexmin, indexmax, tmin, tmax)
    vec = Vector{Vector{Tuple}}()
    j = 10
  
    while true
        if j >= indexmax
            break
        end

        for i in j + indexmin:indexmax
            if matrix[i][j]
                indices = find_orbit(matrix, i, j)
                if length(indices) >= tmin && length(indices) <= tmax && is_original_orbit(vec, indices, indexmax)
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

function find_orbit(matrix, i, j)
    indices = Vector{Tuple}()
    s = size(matrix)[1]

    used = Matrix{Bool}(undef, s, s)
    fill!(used, false)
    
    push!(indices, (i, j))
    used[i, j] = true

    find!(matrix, i, j, indices, used)
    return indices
end

function find!(matrix, i, j, indices, used)
    for i1 in -1:1
        for j1 in -1:1
            if i1 == j1 == 0
                continue
            end
            if matrix[i - i1][j - j1] && !used[i - i1, j - j1]
                push!(indices, (i - i1, j - j1))
                used[i - i1, j - j1] = true
                find!(matrix, i - i1, j - j1, indices, used)
            end
        end
    end
    return nothing
end

function next_column(indices)
    max = 0
    for i in eachindex(indices)
        if (indices[i][2] > max)
            max = indices[i][2]
        end
    end
    return max
end

function is_original_orbit(vec, indices, indexmax)
    n = length(vec)
    if n == 0
        return true
    end

    max = 0
    for i in eachindex(vec[n])
        if (vec[n][i][2] > max)
            max = vec[n][i][2]
        end
    end

    min = indexmax
    for i in eachindex(indices)
        if (indices[i][2] < min)
            min = indices[i][2]
        end
    end

    if max < min
        return true
    end
    return false
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
Y, t = trajectory(lorenz, total_time; Ttr = 3)

R = RecurrenceMatrix(Y, 3.0)

r = Vector{Vector}()
s = size(R)[1]
for i in 0:s:s*(s - 1)
   tmp = Vector{}()
   for j in (i + 1):s + i
       push!(tmp, R[j])
   end
   push!(r, tmp)
end

v = orbits(r, 100, (s - 11), 30, 50)

xs1 = Vector{Int64}()
ys1 = Vector{Int64}()
for i in eachindex(v)
   for j in eachindex(v[i])
       push!(xs1, v[i][j][2])
       push!(ys1, v[i][j][1])
   end
end

CairoMakie.activate!()

current_dir = pwd()

recurrenceplot(R; ascii = true)

fig = Figure(resolution = (10000, 10000))
ax = Axis(fig[1,1])
xs, ys = coordinates(R)
CairoMakie.scatter!(ax, xs, ys; color = :black, markersize = 10)
CairoMakie.scatter!(ax, xs1, ys1; color = :red, markersize = 10)
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

fig2 = Figure(resolution = (1000, 1000))
ax1 = Axis3(fig2[1, 1], azimuth = 0.5 * pi, title = "Original attractor")
GLMakie.lines!(ax1, Y[:, 1], Y[:, 2], Y[:, 3]; linewidth=1, color=:blue)
ax2 = Axis3(fig2[1, 2], azimuth = 0.5 * pi, title = "Reconstructed attractor")
GLMakie.lines!(ax2, x, y, z; linewidth=1, color=:red)
save_path = joinpath(current_dir, "Attractor.png")
save(save_path, fig2)

png_path = joinpath(current_dir, "Png")
rm(png_path, recursive=true, force=true)
mkdir(png_path)

for i in eachindex(v)
    x1 = Float64[]
    y1 = Float64[]
    z1 = Float64[]
    for j in eachindex(v[i])
        push!(x1, Y[v[i][j][1]][1])
        push!(y1, Y[v[i][j][1]][2])
        push!(z1, Y[v[i][j][1]][3])
    end

    figi = Figure(resolution = (1000, 1000))
    axi = Axis3(figi[1, 1], title = "$i orbit")
    GLMakie.lines!(axi, x1, y1, z1; linewidth=5, color=:blue)
    save_path = joinpath(current_dir, "Png\\ReconsructedAttractorOrbit$i.png")
    save(save_path, figi)
end