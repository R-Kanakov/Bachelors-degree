using DynamicalSystems
using CairoMakie
CairoMakie.activate!()

function orbits(matrix, indexmin, indexmax, tmin, tmax)
    s = size(matrix)[1]
    if indexmax > s
        indexmax = s
    end

    vec = Vector{Vector{Tuple}}()
    j = 10

    while true
        if j >= s
            break
        end

        for i in j + indexmin:indexmax 
            if matrix[i][j]
                indices = FindOrbit(matrix, i, j)
                if length(indices) >= tmin && length(indices) <= tmax
                    push!(vec, indices)
                    j = nextcolumn(indices)
                    break
                end
            end
        end

        j += 1
    end
    return vec
end

function find_orbit(matrix, i, j)
    # Создаем вектор пар значений, в котором будем хранить индексы элементов, принадлежащих одной компоненте связности
    indices = Vector{Tuple}()
    s = size(matrix)[1]

    # Инициализируем матрицу, в которой будем помечать какие элементы уже проверялись
    used = Matrix{Bool}(undef, s, s)
    fill!(used, false)
    
    push!(indices, (i, j))
    used[i, j] = true
    # Ищем компоненту связности
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
    v = Vector{}()
    for i in eachindex(indices)
       push!(v, indices[i][2])
    end
    return maximum(v)
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
# Как ни странно R представляет собой не матрицу, а вектор из 4004001 элемента (количество элементов зависит от total_time), поэтому для удобства её необходимо привести в матричный вид
# x[i] - x[j] < 3.0
r = Vector{Vector}()
s = size(R)[1]
for i in 0:s:s*(s - 1)
   tmp = Vector{}()
   for j in (i + 1):s + i
       push!(tmp, R[j])
   end
   push!(r, tmp)
end

v = orbits(r, 100, 1800, 15, 180)

xs1 = Vector{Int64}()
ys1 = Vector{Int64}()
for i in eachindex(v)
   for j in eachindex(v[i])
       push!(xs1, v[i][j][2])
       push!(ys1, v[i][j][1])
   end
end

recurrenceplot(R; ascii = true)
fig = Figure(resolution = (10000, 10000))
ax = Axis(fig[1,1])
xs, ys = coordinates(R)
scatter!(ax, xs, ys; color = :black, markersize = 10)
scatter!(ax, xs1, ys1; color = :red, markersize = 10)
ax.limits = ((1, s), (1, s));
ax.aspect = 1
save("Output.png", fig)