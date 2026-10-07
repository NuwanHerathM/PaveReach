function insert_jacobian(x)
    jacobian_expr = :(Df_fun = ColJacobianMatrix([(x...) -> [(-2*(2.5*sin(x[3]) - x[1]))], (x...) -> [(-2*(2.5*cos(x[3]) - x[2]))], (x...) -> [(5*cos(x[3])*(2.5*sin(x[3]) - x[1]) - 5*sin(x[3])*(2.5*cos(x[3]) - x[2]))]]))
    if x == :({{jacobian}})
        return jacobian_expr
    else
        return x
    end
end

filename = splitext(@__FILE__)[1]

@assert occursin(r"_col$", filename) "The file name does not end with \"_col\"."
fileprefix = filename[begin:(end-4)]
templatefile = "$(fileprefix)_template.jl"

include(insert_jacobian, templatefile)