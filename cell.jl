using AbstractTrees

abstract type Cell end

AbstractTrees.children(cell::Cell) = cell.children
AbstractTrees.ParentLinks(::Cell) = AbstractTrees.StoredParents()
AbstractTrees.parent(cell::Cell) = cell.parent
AbstractTrees.nodevalue(cell::Cell) = cell.box
isleaf(cell::Cell) = isempty(cell.children)
AbstractTrees.isroot(cell::Cell) = isnothing(cell.parent)

function height(cell::Cell)
    if isleaf(cell)
        return 0
    else
        return 1 + maximum(height(child) for child in cell.children)
    end
end

#-----------------------------------------------------------------------------------------

mutable struct ApproximationCell <: Cell
    const box::IntervalBox
    children::Vector{ApproximationCell}
    parent::Union{ApproximationCell, Nothing}
    inner
    outer
end

function Base.show(io::IO, cell::ApproximationCell)
    println("Inner: [", cell.inner.lo, ", ", cell.inner.hi, "]")
    print("Outer: [", cell.outer.lo, ", ", cell.outer.hi, "]")
end

function tolerance(cell::ApproximationCell)
    return max(abs(cell.inner.lo - cell.outer.lo), abs(cell.inner.hi - cell.outer.hi))
end

function make_approximationcell_root(box::IntervalBox, inner, outer)
    return ApproximationCell(box, [], nothing, inner, outer)
end

function make_approximationcell_leaf(box::IntervalBox, inner, outer, parent::Cell)
    cell = ApproximationCell(box, [], parent, inner, outer)
    push!(parent.children, cell)
    return cell
end

#-----------------------------------------------------------------------------------------

mutable struct MembershipCell <: Cell
    const box::IntervalBox
    children::Vector{MembershipCell}
    parent::Union{MembershipCell, Nothing}
    is_in::Function
    is_out:: Function
end

function make_membershipcell_root(box::IntervalBox, is_in, is_out)
    return MembershipCell(box, [], nothing, is_in, is_out)
end

function make_membershipcell_leaf(box::IntervalBox, is_in, is_out, parent::Cell)
    cell = MembershipCell(box, [], parent, is_in, is_out)
    push!(parent.children, cell)
    return cell
end

function diam(cell::MembershipCell)
    if isleaf(cell)
        return IntervalArithmetic.diam(cell.box)
    else
        return maximum(diam.(cell.children))
    end
end
