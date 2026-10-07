struct UserDefinedFunctions
    f::Vector{Function}
    Df::JacobianMatrix
end

get_f(userfunctions::UserDefinedFunctions) = userfunctions.f
get_Df(userfunctions::UserDefinedFunctions) = userfunctions.Df

struct UserDefinedSymbolicFunctions
    variables::AbstractVector{Num}
    f_num::Vector{Num}
end

get_variables(userfunctions::UserDefinedSymbolicFunctions) = userfunctions.variables
get_f_num(userfunctions::UserDefinedSymbolicFunctions) = userfunctions.f_num
