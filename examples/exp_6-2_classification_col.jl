function insert_jacobian(x)
    jacobian_expr = :(Df_fun = ColJacobianMatrix(collect((x...) -> [confidence(dN(i, x[1:p]))] for i in 1:p)))
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