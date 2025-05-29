#using Pkg; Pkg.add(url = "https://github.com/JuliaDynamics/PeriodicOrbits.jl")
#using Optim
#using PeriodicOrbits
#using Base.Threads
#using IntervalTrees

# TODO: try CairoMakie.activate!(type = "svg")

using DynamicalSystems
using DifferentialEquations
using CairoMakie
using GLMakie
using LinearAlgebra
using ArgParse
using SparseArrays
using Statistics

global s
global png_path
global p0

function sparse_dilation(cols, rows)
    new_rows = Int[]
    new_cols = Int[]
    
    offsets = [-1  0;  # Up
                1  0;  # Down
                0 -1;  # Left
                0  1;  # Right
               -1 -1;  # Up-left
               -1  1;  # Up-right
                1 -1;  # Down-left
                1  1]  # Down-right
    
    for (r, c) in zip(rows, cols)
        for offset in eachrow(offsets)
            new_r = r + offset[1]
            new_c = c + offset[2]
            
            if  1 <= new_r <= s && 1 <= new_c <= s
                push!(new_rows, new_r)
                push!(new_cols, new_c)
            end
        end
        
        push!(new_rows, r)
        push!(new_cols, c)
    end
    
    unique_indices = unique(zip(new_rows, new_cols))
    
    final_rows = [i[1] for i in unique_indices]
    final_cols = [i[2] for i in unique_indices]
    
    return final_cols, final_rows
end

"""
Optimization
From 
total_time = 500   -> 180 s
To
total_time = 1000  -> 0.38 s
total_time = 10000 -> 0.80 s
"""
function orbits(R, t_max, t_min, dt::Float64; save_matrix = save_matrix)
    if p == 0
        x, y = coordinates(R)
        save_matrix("1_Base", x, y, (ax)->())
    end

    # Delete lower part of matrix
    R_sparce = sparse(RecurrenceAnalysis.tril(R, -1))

    R_rows = rowvals(R_sparce)
    R_cols = findnz(R_sparce)[2]
    
    p == 0 && save_matrix("2_Without_Bottom", R_cols, R_rows, (ax)->())

    # Removing all points that don't satisfy the
    # conditions of being between t_min and t_max
    times = [(R_rows[i] - R_cols[i]) * dt for i in eachindex(R_cols)]

    mask = t_min .<= times .<= t_max

    if p == 0
        function add_lines(ax)
            CairoMakie.lines!(ax, [0, s - t_min / dt], [t_min / dt, s];
                              color = :red)
            CairoMakie.lines!(ax, [0, s - t_max / dt], [t_max / dt, s];
                              color = :red)
        end
        save_matrix("3_Time_Limits", R_cols, R_rows, add_lines)
    end

    R_rows = R_rows[mask]
    R_cols = R_cols[mask]

    p == 0 && save_matrix("4_Time_Limits_Del", R_cols, R_rows, add_lines)

    # Dilating the matrix
    R_cols, R_rows = sparse_dilation(R_cols, R_rows)

    p == 0 && save_matrix("5_With_Dilation", R_cols, R_rows, (ax)->())

    # Removing all points that are higher than one of the found point
    R_cols_len = length(R_cols)

    pairs = [(R_rows[i], R_cols[i]) for i in eachindex(R_rows)]

    sort!(pairs, by = x -> (x[2], x[1]))

    unique_pairs = Set{}()
    mask = falses(length(pairs))

    for (i, pair) in enumerate(pairs)
        if pair[2] ∉ unique_pairs
            push!(unique_pairs, pair[2])
            mask[i] = true
        end
    end

    R_rows = [pairs[i][1] for i in eachindex(pairs) if mask[i]]
    R_cols = [pairs[i][2] for i in eachindex(pairs) if mask[i]]

    p == 0 && save_matrix("6_Without_Duplicates", R_cols, R_rows, (ax)->())

    # Removing all points that are to the right of one of the found point
    R_cols_len = length(R_cols)
    mask = trues(R_cols_len)

    for i in 1:R_cols_len - 1
        println("$i, R_cols[i] = $(R_cols[i]), R_cols[i + 1] = $(R_cols[i + 1]), R_rows[i] = $(R_rows[i]), R_rows[i + 1] = $(R_rows[i + 1])")
        if R_cols[i] + 1 == R_cols[i + 1] &&
           R_rows[i] + 1 == R_rows[i + 1]
            mask[i + 1] = false
        elseif R_cols[i] + 1 == R_cols[i + 1] &&
               R_rows[i] == R_rows[i + 1]
            mask[i + 1] = false
        end
    end

    R_rows = R_rows[mask]
    R_cols = R_cols[mask]

    p == 0 && save_matrix("7_Only_Dots", R_cols, R_rows, (ax)->())

    length(R_rows) == 0 && (error("No trajectories lying so close have been found.
                                   Try to increase total_time or ε"))

    return R_rows, R_cols
end

@inbounds function lorenz_rule!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
    return nothing
end

function run_simulation(choice::Int64, total_time::Float64,
                        u0::Vector{Float64},
                        ε::Float64, t_max, t_min)
    if choice == 0
        dt = 0.01

        lorenz = CoupledODEs(lorenz_rule!, u0, p0)

        Y, t = trajectory(lorenz, total_time; Ttr = 3, Δt = dt)
    elseif choice == 1
        tspan = (0.0, total_time)

        prob = ODEProblem(lorenz_rule!, u0, tspan, p0)
        sol = solve(prob, Tsit5())

        Y = StateSpaceSet(sol.u)
        t = sol.t

        dt = round(mean(diff(t)), digits=3)
    end

    R = RecurrenceMatrix(Y, ε)

    global s = length(t)

    R_rows, R_cols = orbits(R, t_max, t_min, dt)

    return Y, R, R_rows, R_cols
end

function save_matrix(name, cols, rows, f)
    # Recurrence matrix visualization
    CairoMakie.activate!()
    
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], title = "Recurrence matrix")
    CairoMakie.scatter!(ax, cols, rows; color = :black)
    f(ax)
    ax.limits = ((1, s), (1, s))
    ax.aspect = 1
    save_path = joinpath(png_path, "Recurrence_Matrix_$name.png")
    save(save_path, fig)
end

function save_matrix_(name, cols, rows, f)
    # Recurrence matrix visualization
    CairoMakie.activate!()
    
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], title = "Recurrence matrix", titlesize = 35, xticklabelsize = 30, yticklabelsize = 30)
    CairoMakie.scatter!(ax, cols, rows; color = :black)
    f(ax)
    ax.limits = ((1, 100), (1, 100))
    ax.aspect = 1
    save_path = joinpath(png_path, "Recurrence_Matrix_$name.png")
    save(save_path, fig)
end


function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--choice", "-c"
            help = "0 to use a continuous dynamic system solution or 1 to
                    use an adaptive step solution"
            arg_type = Int
            default = 0
        "--total_time", "-t"
            help = "Total time of dynamic system evaluation"
            arg_type = Float64
            default = 50
        "--ε", "-e"
            help = "Parameter for recurrence matrix"
            arg_type = Float64
            default = 0.5
        "--t_max"
            help = "Maximum time of first return"
            default = 10
        "--t_min"
            help = "Minimum time of first return"
            default = 0.05
        "--p", "-p"
            help = "0 to plot all steps of matrix transform
                    1 to skip this"
            arg_type = Int64
            default = 0
    end

    parsed_args = parse_args(s)

    choice = parsed_args["choice"]
    total_time = parsed_args["total_time"]
    ε = parsed_args["ε"]
    t_max = parsed_args["t_max"]
    t_min = parsed_args["t_min"]
    global p = parsed_args["p"]

    if choice ∉ [0, 1]
        error("Invalid choice. Please select 1 to use a continuous dynamic 
               system solution or 2 to use an adaptive step solution")
    end

    if total_time <= 0
        error("Must be greater then 0.")
    end

    return choice, total_time, ε, t_max, t_min
end

function main()
    # Making directories and paths to them
    current_dir = dirname(dirname(dirname(@__FILE__)))
    # For png
    global png_path = joinpath(current_dir, "Png","Png_old")
    rm(png_path, recursive = true, force = true)
    mkdir(png_path)
    # For mkv
    mkv_path = joinpath(current_dir, "Mkv", "Mkv_old")
    rm(mkv_path, recursive = true, force = true)
    mkdir(mkv_path)

    choice, total_time, ε, t_max, t_min = parse_arguments()

    u0 = [10, 10, 10.9]
    global p0 = [10, 28, 2.666]

    #Y, R, R_rows, R_cols = run_simulation(choice, total_time, u0, ε, t_max, t_min)

    n = 100
    global s = 100
    R_data = zeros(n, n)

    for i in 1:n
        for j in 1:n
            if i == j || abs(i - j) <= 5 || (i % 10 == 0 && j % 10 == 0)
                R_data[i, j] = 1
            end
        end
    end
    R_data[21,11]=1
    
    orbits(sparse(R_data), 0.5, 0.09, 0.01; save_matrix = save_matrix_)

    @assert 1 == 2

    # Visualization and animation
    x, y = coordinates(R)
    scat = (ax) -> CairoMakie.scatter!(ax, R_cols, R_rows; color = :red, markersize = 12)
    save_matrix("0", x, y, scat)
    
    GLMakie.activate!()

    # Array which will represent distance between
    # start and last dots of periodic orbit in indexes
    indexDistance = [[R_cols[i], R_rows[i]] for i in eachindex(R_rows)]

    # Attractor visialization and orbits animation
    xmin, xmax, ymin, ymax, zmin, zmax =
            minimum(p[1] for p in Y), maximum(p[1] for p in Y),
            minimum(p[2] for p in Y), maximum(p[2] for p in Y),
            minimum(p[3] for p in Y), maximum(p[3] for p in Y)

    fig2 = Figure(size=(5000, 2500))
    ax2 = Axis3(fig2[1, 1], azimuth = 0.5 * pi, titlesize = 40, title = "Original attractor")
    limits!(ax2, xmin, xmax, ymin, ymax, zmin, zmax)
    ax3 = Axis3(fig2[1, 2], azimuth = 0.5 * pi, titlesize = 40, title = "Reconstructed attractor")
    limits!(ax3, xmin, xmax, ymin, ymax, zmin, zmax)
    GLMakie.lines!(ax2, Y[:, 1], Y[:, 2], Y[:, 3]; linewidth = 1, color = :blue)

    colors = [:blue, :red, :green, :orange, :black, :purple, :pink,
              :brown, :cyan, :magenta, :turquoise, :coral, :lavender,
              :salmon, :teal, :gray, :violet, :gold, :silver, :yellow]

    k, n = 0, length(colors)

    for i in eachindex(indexDistance)

        points = [Y[j] for j in indexDistance[i][1]:indexDistance[i][2]]

        # TODO: Check orbit with ? method

        # Reconstructed attractor visialization
        k = (i - 1) % n + 1
        GLMakie.lines!(ax3, points; linewidth = 1, color = colors[k])

        # Orbits animation
        #fig3 = Figure(size=(1000, 1000))
        #ax4 = Axis3(fig3[1, 1])
        #limits!(ax4, xmin, xmax, ymin, ymax, zmin, zmax)
        # save_path = joinpath(mkv_path, "ReconsructedAttractorOrbit$i.mp4")
        # println(length(points))
        # @time begin
        #     record(fig3, save_path) do io
        #         pointObservable = Observable(Point3f[points[1]])
        #         for j in 2:length(points)
        #             pointObservable[] = push!(pointObservable[], points[j])
        #             GLMakie.lines!(ax4, pointObservable[]; color = :blue)#, markersize = 5000)
        #             recordframe!(io)
        #         end
        #     end
        # end
    end

    save_path = joinpath(png_path, "Attractor.png")
    save(save_path, fig2)
    return nothing
end

main()
