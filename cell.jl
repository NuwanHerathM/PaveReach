using AbstractTrees

abstract type AbstractCell end

struct Cell <: AbstractCell
    interval::IntervalArithmetic.Interval{T} where T<:Number
    children::Vector{AbstractCell}
    parent::AbstractCell
end

struct CellStart <: AbstractCell
    children::Vector{AbstractCell}
end

struct CellEnd <: AbstractCell
    interval::IntervalArithmetic.Interval{T} where T<:Number
    parent::AbstractCell
    is_in::Function
    is_out::Function
end

AbstractTrees.children(cell::AbstractCell) = isa(cell, CellEnd) ? [] : cell.children
AbstractTrees.parent(cell::AbstractCell) = isa(cell, CellStart) ? nothing : cell.parent
AbstractTrees.nodevalue(cell::AbstractCell) = isa(cell, CellStart) ? emptyinterval() : cell.interval
isleaf(cell::AbstractCell) = isa(cell, CellEnd) || (isroot(cell) && isempty(cell.children))
AbstractTrees.isroot(cell::AbstractCell) = isa(cell, CellStart)

function Cell(interval::IntervalArithmetic.Interval{T}, parent::Union{CellStart, Cell}) where T<:Number
    return Cell(interval, [], parent)
end

function CellStart()
    return CellStart([])
end

import Base: push!

function push!(cell::CellStart, intervals, is_in, is_out)
    return push!(cell, intervals, 1, is_in, is_out)
end

function push!(cell::AbstractCell, intervals, i, is_in, is_out)
    if i > length(intervals)
        return
    end

    pos = findfirst(isequal(intervals[i]), [child.interval for child in cell.children])
    if !isnothing(pos)
        push!(cell.children[pos], intervals, i+1, is_in, is_out)
    else
        child = make_paving(intervals, i, cell, is_in, is_out)
        push!(cell.children, child)
    end
end

function remove!(cell::AbstractCell)
    if isroot(cell)
        return
    end

    parent = cell.parent
    filter!(child -> child !== cell, parent.children)

    if isempty(parent.children)
        remove!(parent)
    end
end

function remove!(cell::CellEnd, intervals)
    if isempty(intervals)
        remove!(cell)
    end
end

function remove!(cell::Union{CellStart, Cell}, intervals)
    if isempty(intervals)
        return
    end
    
    pos = findfirst(isequal(intervals[begin]), [child.interval for child in cell.children])
    if !isnothing(pos)
        remove!(cell.children[pos], intervals[2:end])
    end
end

function diam(cell::AbstractCell)
    if isleaf(cell)
        return IntervalArithmetic.diam(cell.interval)
    elseif isroot(cell)
        return maximum(diam.(cell.children))
    else
        return maximum([IntervalArithmetic.diam(cell.interval), diam.(cell.children)...])
    end
end

#-----------------------------------------------------------------------------------------

mutable struct ApproximationCell
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


