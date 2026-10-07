function insert_jacobian(x)
    # General formula
    # df_1(x) = [confidence_1(dN(x[1:2], 1)), confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 1]), confidence_1(dN(x[1:2], 1)), confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 1])]
    # df_2(x) = [confidence_1(dN(x[1:2], 2)), confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 2]), confidence_1(dN(x[1:2], 2)), confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 2])]
    # df_3(x) = [0, confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 3]), 0, confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 3])]
    # df_4(x) = [0, confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 4]), 0, confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 4])]
    # Formula in practice
    df_1(x) = [confidence_1(dN(x[1:2], 1)), confidence_1(dN(perturbation(x[1:4]), 1)), confidence_1(dN(x[1:2], 1)), confidence_1(dN(perturbation(x[1:4]), 1))]
    df_2(x) = [confidence_1(dN(x[1:2], 2)), confidence_1(dN(perturbation(x[1:4]), 2)), confidence_1(dN(x[1:2], 2)), confidence_1(dN(perturbation(x[1:4]), 2))]
    df_3(x) = [0, confidence_1(dN(perturbation(x[1:4]), 1)), 0, confidence_1(dN(perturbation(x[1:4]), 1))]
    df_4(x) = [0, confidence_1(dN(perturbation(x[1:4]), 2)), 0, confidence_1(dN(perturbation(x[1:4]), 2))]
    jacobian_expr = :(Df_fun = ColJacobianMatrix([df_1, df_2, df_3, df_4]))
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