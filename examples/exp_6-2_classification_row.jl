function insert_jacobian(x)
    jacobian_expr = :(Df_fun = RowJacobianMatrix([x -> [confidence(DN(IntervalBox([x[1], x[2], x[3]])))...]]))
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