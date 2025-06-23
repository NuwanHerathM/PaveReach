include("genreach.jl")

using Statistics

@variables x[1:3]
QE5=([x[1]^2/4.0+(x[2]+1.0)*(x[3]+2.0)+(x[3]+3.0)^2],["exists", 1, "forall", 2, "exists", 3],3,1)

g = QE5[1]

range = -1:0.01:1
l = length(range)
out = Array{Union{Missing, Float64}}(missing, l, l, l)
p = 3
n = 1

for j in 1:n
    g_expr = build_function(g[j], [x[i] for i=1:p])
    my_g = eval(g_expr)
    idx_x = 1
    for i in range
        idx_y = 1
        for j in range
            idx_z = 1
            for k in range
                input_center = [i, j, k]
                out[idx_x, idx_y, idx_z] = Base.invokelatest(my_g,input_center)
                idx_z += 1
            end
            idx_y += 1
        end
        idx_x += 1
    end
end

idx_z_min = argmin(out, dims = 3)
idx_z_max = argmax(out, dims = 3)
idx_y_min = argmax(out[idx_z_min], dims = 2)
idx_y_max = argmin(out[idx_z_max], dims = 2)
println("Min: ", minimum(out[idx_z_min][idx_y_min]))
println("Max: ", maximum(out[idx_z_max][idx_y_max]))