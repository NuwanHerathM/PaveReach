function insert_jacobian(x)
    jacobian_expr = :(Df_fun = RowJacobianMatrix([(x...) -> [(-2*(2.5*sin(x[3]) - x[1])), (-2*(2.5*cos(x[3]) - x[2])), (5*cos(x[3])*(2.5*sin(x[3]) - x[1]) - 5*sin(x[3])*(2.5*cos(x[3]) - x[2]))]]))
    if x == :({{jacobian}})
        return jacobian_expr
    else
        return x
    end
end

filename = splitext(@__FILE__)[1]

@assert occursin(r"_row$", filename) "The file name does not end with \"_row\"."
fileprefix = filename[begin:(end-4)]
templatefile = "$(fileprefix)_template.jl"

include(insert_jacobian, templatefile)