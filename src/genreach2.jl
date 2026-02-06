using IntervalArithmetic
using Symbolics 
using LazySets

function print_QE(QE)
  print("Computing approximations for R={z | ")

 for k in 1:2:length(QE[2])-1
     print(QE[2][k])
     print(" x[")
     print(QE[2][k+1])
     print("] ")
 end
 print("z = (")
 for j in 1:QE[4]-1
   print(QE[1][j])
   print(", ")
 end
 print(QE[1][QE[4]])
 println(")}")
end

function _isapprox_hack(x,y)
  return isapprox(x,y)
end

function isapproxzero_hack(x)
  return x==0.0
end

function ismultiple_hack(u::AbstractVector{<:Real}, v::AbstractVector{<:Real})
  return _ismultiple_hack(u, v; allow_negative=true)
end

function _ismultiple_hack(u::AbstractVector, v::AbstractVector; allow_negative::Bool)
  @assert length(u) == length(v) "wrong dimension"
  no_factor = true
  factor = 0
  @inbounds for i in eachindex(u)
      if isapproxzero_hack(u[i])
          if !isapproxzero_hack(v[i])
              return (false, 0)
          end
          continue
      elseif isapproxzero_hack(v[i])
          return (false, 0)
      end
      if no_factor
          no_factor = false
          factor = u[i] / v[i]
          if !allow_negative && factor < 0
              return (false, 0)
          end
      elseif !_isapprox_hack(factor, u[i] / v[i])
          return (false, 0)
      end
  end
  if no_factor
      # both vectors are zero
      return (true, 0)
  end
  return (true, factor)
end

function nonzero_columns_hack(A::AbstractMatrix; comparison=isapproxzero_hack)
  n = size(A, 2)
  nzcol = Vector{Int}()
  sizehint!(nzcol, n)
  for j in 1:n
      if !all(comparison, view(A, :, j))
          push!(nzcol, j)
      end
  end
  return nzcol
end

function remove_zero_columns_hack(A::AbstractMatrix)
  nzcol = nonzero_columns_hack(A)
  if length(nzcol) == size(A, 2)
      return A
  else
      return A[:, nzcol]
  end
end

function remove_redundant_generators_hack(Z::Zonotope{N}) where {N}
  G = genmat(Z)
  G = remove_zero_columns_hack(G)
  p = size(G, 2)
  removed_zero_generators = p < ngens(Z)
  deleted = false
  done = falses(p)
  #G_new = _vector_type(typeof(G))[]  # list of new column vectors
  G_new = Vector{N}[]
  @inbounds for j1 in 1:p
      if done[j1]  # skip if the generator was already removed
          continue
      end
      # "done[j1] = true" not needed because we will never look at it again
      gj1 = G[:, j1]
      for j2 in (j1+1):p  # look at all generators to the right
          if done[j2]  # skip if the generator was already removed
              continue
          end
          gj2 = G[:, j2]
          answer, factor = ismultiple_hack(gj1, gj2)
          if answer
              # column j2 is a multiple of column j1
              if factor > zero(N)
                  gj1 += gj2
              else
                  gj1 -= gj2
              end
              done[j2] = true
              deleted = true
          end
      end
      push!(G_new, gj1)
  end

  if deleted
      G_new = reduce(hcat, G_new)  # convert list of column vectors to matrix
      return Zonotope(center(Z), G_new)
  elseif removed_zero_generators
      return Zonotope(center(Z), G)
  end
  return Z  # return the original zonotope if no generator was removed
end

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

function print_vol(QE)
  try
    global voli=LazySets.volume(QE[1])
  catch e
    global voli=0
  end
  try 
    global volo=LazySets.volume(QE[2])
  catch e
    global volo=0
  end
  print("Volume inner=")
  println(voli)
  print("Volume outer=")
  println(volo)
  print("Relative precision=")
  println(voli/volo)
end

function print_vol_o0(res, p)
  global voli = 1.0
  global volo = 1.0

  for i in 1:p
    try 
      global voli = voli*diameter(res[1][i])
    catch e
      global voli = 0
    end
    try
      global volo = volo*diameter(res[2][i])
    catch e
      global volo = 0
    end
    
  end
  
  print("Volume inner=")
  println(voli)
  print("Volume outer=")
  println(volo)
  print("Relative precision=")
  println(voli/volo)
end

function print_rel_vol_o1_o0(res1, res0, p)
  try 
    global voli1=LazySets.volume(res1[1])
  catch e
    global voli1=0
  end
  try
    global volo1=LazySets.volume(res1[2])
  catch e
    global volo1=0
  end
  
  
  global voli0 = 1.0
  global volo0 = 1.0

  for i in 1:p
    try
      global voli0 = voli0*diameter(res0[1][i])
    catch e
      global voli0 = 0
    end
    try
      global volo0 = volo0*diameter(res0[2][i])
    catch e
      global volo0 = 0
    end
    
  end

  print("Vol(outer1)/Vol(outer0)=")
  println(volo1/volo0)
  println

  print("Vol(inner1)/Vol(inner0)=")
  println(voli1/voli0)
  println
end

function Zon(lin, k, n)
  no=0.0
  for j in 1:n
    no = no+abs(lin[j,k])
  end
  if no<0.0000001 # hack horrible pour eviter plantage de minkowski_sum sur les polytopes
    return ZeroSet(n)
  else
  return Zonotope(zeros(n),reshape([lin[j,k] for j=1:n],n,1))
  end
end

function IZon(lin, k, n)
  return LazySets.Interval(-abs(lin[1,k]),abs(lin[1,k]))   
end

function build_function_f_Df(f_num::Vector{Num}, x, n::Int, p::Int)
  f_fun = Function[]
  Df_fun = Function[]
  for j in 1:n
    f_expr = build_function(f_num[j], [x[i] for i=1:p])
    push!(f_fun, eval(f_expr))
    
    Df = Symbolics.jacobian([f_num[j]], [x[i] for i=1:p])
    Df_expr = build_function(Df, [x[i] for i=1:p])
    push!(Df_fun, eval(Df_expr[1]))
  end
  return f_fun, Df_fun
end

function QEapprox_o0(g_fun, Dg_fun, quantifiers, q, p, n, input)
  global out = [] #LazySets.Interval(0.0,0.0) for k = 1:n] 
  global inn = [] #LazySets.Interval(0.0,0.0) for k = 1:n] 

  radii = [IntervalArithmetic.radius(input[i]) for i=1:p]
  input_center = [mid(input[i]) for i=1:p]
   
  for j in 1:n
     c = Base.invokelatest(g_fun[j],input_center)
     range_Dg = Base.invokelatest(Dg_fun[j],input)

     outer = interval(c,c)
      # println(outer)
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

    inner = interval(c,c)
  
    for i = 2p-1:-2:1
      if q[j][i] == "exists"
        inner = inner+I(range_Dg,radii,q[j][i+1])
      else
        inner = slash(inner,O(range_Dg,radii,q[j][i+1]))
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
  
  return (inn, out)
end

function QEapprox_o1(g, quantifiers, p, n)
input = [interval(-1.0,1.0) for i = 1:p]   # @SVector but pb with scope 
input_center = zeros(p)
input_t=transpose(input)
  
linear = Array{Float64,2}(undef,n,p)
range_H = Array{IntervalArithmetic.Interval,1}(undef,n) # Array{IntervalArithmetic.Interval,1}(undef,n)
c = Array{Float64,1}(undef,n)

for j in 1:n
  g_expr = build_function(g[j], [x[i] for i=1:p])
  my_g = eval(g_expr)
  c[j] = Base.invokelatest(my_g,input_center)
  
  # Jacobian in 1xp matrix form
  Dg = Symbolics.jacobian([g[j]], [x[i] for i=1:p])
  Dg_expr = build_function(Dg, [x[i] for i=1:p])
  my_Dg = eval(Dg_expr[1])
  for k in 1:p
    linear[j,k] = Base.invokelatest(my_Dg,input_center)[k]
  end
  
  # Hessian evaluation on the full hypercube
  Hg = Symbolics.hessian(g[j], [x[i] for i=1:p]) # was [g[j]]
  H_expr = build_function(Hg, [x[i] for i=1:p])
  my_H = eval(H_expr[1])
  range_Hg = Base.invokelatest(my_H,input)
  range_H[j] = 0.0
  for i in 1:p
    range_H[j] = range_H[j]+range_Hg[i,i]*interval(0.0,1.0)
    for k in (i+1):p
      range_H[j] = range_H[j]+2*range_Hg[i,k]*interval(-1.0,1.0)
    end
  end
  range_H[j] = 0.5*range_H[j]
end
# Calcul dans le cas linéaire (sur linear[j]+u)
  # outerapprox first
  u = [diam(range_H[i])/2.0 for i=1:n] 
  global uez = true
  global outer
  global inner
  for i in 1:n
    if u[i]!==0.0
      global uez = false
    end
  end  
  
  # particular case for n=1
  if n==1
    Icu = mid(range_H[1])
    if uez
      global outer = ZeroSet(1)
    else
      global outer = LazySets.Interval(-u[1],u[1])
    end
  
    for i = 2p-1:-2:1
      IZo = IZon(linear,quantifiers[i+1],n)
      if quantifiers[i] == "exists"
         global outer = minkowski_sum(outer,IZo)
      else
         global outer = minkowski_difference(outer,IZo)
         if (outer==EmptySet(p))
          break
         end
      end
    end

    global outer = LazySets.translate(outer,[c[1].+Icu])

  # inner approximation
  
    global inner = ZeroSet(1)

    for i = 2p-1:-2:1
      IZo = IZon(linear,quantifiers[i+1],n)
      if quantifiers[i] == "exists"
         global inner = minkowski_sum(inner,IZo) 
      else
         # emptyness condition to be done...
         global inner = minkowski_difference(inner,IZo)
         if (inner==EmptySet(p))
          break
         end
      end
    end

    # last quantifier (forall, for the u part)
    if !uez
      global inner = minkowski_difference(inner,LazySets.Interval(-u[1],u[1]))
    end
    
    if (inner!=EmptySet(p))
      global inner = LazySets.translate(inner,[c[1].+Icu])  
    end

    return (inner, outer)
  else
  # general case n>=1
    cu = [mid(range_H[i]) for i=1:n]
    if uez
      global outer = ZeroSet(n)
    else
      global outer = Zonotope(zeros(n),reshape(u,n,1))
    end
  
    for i = 2p-1:-2:1
      Zo = Zon(linear,quantifiers[i+1],n)
      if quantifiers[i] == "exists"
         global outer = minkowski_sum(outer,Zo,prune=false)
      else
         try
          global outer = remove_redundant_generators_hack(outer)
         catch e
         end
         global outer = minkowski_difference(outer,Zo)
         if (outer==EmptySet(p))
          break
         end
      end
    end

    global outer = LazySets.translate(outer,c.+cu)

  # inner approximation
  
    global inner = ZeroSet(n)

    for i = 2p-1:-2:1
      Zo = Zon(linear,quantifiers[i+1],n)
      if quantifiers[i] == "exists"
         global inner = minkowski_sum(inner,Zo,prune=false) 
      else
         try
          global inner = remove_redundant_generators_hack(inner) # + remove_zero_generators
         catch e
         end
         global inner = minkowski_difference(inner,Zo)
         if (inner==EmptySet(p))
          break
         end
      end
   end

# last quantifier (forall, for the u part)
    if !uez
      try
        global inner = remove_redundant_generators_hack(inner)
      catch e
      end
      global inner = minkowski_difference(inner,Zonotope(zeros(n),reshape(u,n,1)))
    end
    
    if (inner!=EmptySet(p))
      global inner = LazySets.translate(inner,c.+cu)
    end

    return (inner, outer) # [inner outer]
  end  
end 


