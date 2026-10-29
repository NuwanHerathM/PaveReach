# Functions extracted from https://github.com/goubault/GenReach

function slash(I, J)
  # convenience function for interpretation of "forall"
  # empty case taken care of here - otherwise use emptyinterval()
  if I.lo+J.hi<=I.hi+J.lo
    return interval(I.lo+J.hi, I.hi+J.lo)
  else
    return EmptySet(1)
  end
end
  
  function O(range_Jf, radii, i)
  # convenience function returning the outer-approximated interval for the contribution of xi, i between 1 and dim input, to function f for which we give the range of the Jacobian range_Jf
    return (range_Jf[i])*interval(-1, 1)*radii[i]
  end
  
  function I(range_Jf, radii, i)
  # convenience function returning the inner-approximated interval for the contribution of xi, i between 1 and dim input, to function f for which we give the Jacobian Jf
  # try catch because the Jacobian of a linear function gives Float64 and not Interval
    try
      Jl = abs(range_Jf[i]).lo
      return Jl*interval(-1, 1)*radii[i]
      catch u
        # println(u)
        Jl = abs(range_Jf[i])
        res = Jl*interval(-1, 1)
        return res*radii[i]
    end
  end

function build_function_f_Df(f_num, x, n::Int, p::Int)
  f_fun = Function[]
  for i in 1:n
    f_expr = build_function(f_num[i], x...)
    push!(f_fun, eval(f_expr))
  end
  Df_fun = Matrix{Function}(undef, n, p)
  Df_num = Symbolics.jacobian(f_num, x)
  for i in 1:n
    for j in 1:p
      Df_expr = build_function(Df_num[i, j], x...)
      Df_fun[i, j] = eval(Df_expr)
    end
  end
  return f_fun, Df_fun
end

function QEapprox_o0_inner(g_fun, Dg_fun, quantifiers_relaxed, n, input)
  global inn = [] #LazySets.Interval(0.0,0.0) for k = 1:n] 
  
  p = length(input)
  radii = IntervalArithmetic.radius.(input)
  input_center = mid.(input)
  
  inputcs = []
  for k in 1:p
    inputc = Vector{IntervalArithmetic.Interval{Float64}}(undef, p)
    inputc[1:k] = input[1:k]
    inputc[k+1:p] = [interval(c, c) for c in input_center[k+1:p]]
    push!(inputcs, inputc)
  end

  for j in 1:n
    c = Base.invokelatest(g_fun[j], input_center...)
    # modif for the refined generalized mean value theorem
    #    range_Dg = Base.invokelatest(Dg_fun[j],input)
    
    range_Dg = Vector(undef,p)
    for k in 1:p
      # Dg_j = Base.invokelatest(Dg_fun[j], inputc)
      # range_Dg[k] = Dg_j[k]
      range_Dg[k] = Base.invokelatest(Dg_fun[j,k], inputcs[k]...)
    end
    
    inner = interval(c,c)
    
    for i = 2p-1:-2:1
      if quantifiers_relaxed[j][i] == "exists"
        inner = inner+I(range_Dg,radii,quantifiers_relaxed[j][i+1])
      else
        inner = slash(inner,O(range_Dg,radii,quantifiers_relaxed[j][i+1]))
        if (inner==EmptySet(1))
          break
        end
      end
    end
    
    if (inner==EmptySet(1))
      push!(inn,inner) 
    else    
      push!(inn,LazySets.Interval(inner))
    end
  end
  
  return inn
end

function QEapprox_o0_outer(g_fun, Dg_fun, quantifiers, n, input)
  global out = [] #LazySets.Interval(0.0,0.0) for k = 1:n] 
  
  p = length(input)
  radii = IntervalArithmetic.radius.(input)
  input_center = mid.(input)
  
  inputcs = []
  for k in 1:p
    inputc = Vector{IntervalArithmetic.Interval{Float64}}(undef, p)
    inputc[1:k] = input[1:k]
    inputc[k+1:p] = [interval(c, c) for c in input_center[k+1:p]]
    push!(inputcs, inputc)
  end

  for j in 1:n
    c = Base.invokelatest(g_fun[j], input_center...)
    # modif for the refined generalized mean value theorem
    #    range_Dg = Base.invokelatest(Dg_fun[j],input)
    
      range_Dg = Vector(undef,p)
      for k in 1:p
        # Dg_j = Base.invokelatest(Dg_fun[j], inputc)
        # range_Dg[k] = Dg_j[k]
          range_Dg[k] = Base.invokelatest(Dg_fun[j,k], inputcs[k]...)
      end
    
    outer = interval(c,c)

    for i = 2p-1:-2:1
      if quantifiers[i] == "exists"
        outer = outer+O(range_Dg,radii,quantifiers[i+1])
      else
        outer = slash(outer,I(range_Dg,radii,quantifiers[i+1]))
        if (outer==EmptySet(1))
          break
        end
      end
      # println(outer)
    end
  
    if (outer==EmptySet(1))
      push!(out,outer) 
    else
      push!(out,LazySets.Interval(outer))
    end
  end
  
  return out
end
