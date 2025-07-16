using AbstractTrees

abstract type AbstractCell end

abstract type Cell <: AbstractCell end

struct BisectableCell <: Cell
    interval::IntervalArithmetic.Interval{T} where T<:Number
    children::Vector{AbstractCell}
    parent::AbstractCell
end

struct AtomicCell <: Cell
    interval::IntervalArithmetic.Interval{T} where T<:Number
    children::Vector{AbstractCell}
    parent::AbstractCell
end

struct CellStart <: AbstractCell
    children::Vector{AbstractCell}
end

abstract type CellEnd <: AbstractCell end

struct BisectableCellEnd <: CellEnd
    interval::IntervalArithmetic.Interval{T} where T<:Number
    parent::AbstractCell
    is_member::Function
end

struct AtomicCellEnd <: CellEnd
    interval::IntervalArithmetic.Interval{T} where T<:Number
    parent::AbstractCell
    is_member::Function
end

Bisectable = Union{BisectableCell, BisectableCellEnd}
Atomic = Union{AtomicCell, AtomicCellEnd}

AbstractTrees.children(cell::AbstractCell) = isa(cell, CellEnd) ? [] : cell.children
AbstractTrees.parent(cell::AbstractCell) = isa(cell, CellStart) ? nothing : cell.parent
AbstractTrees.nodevalue(cell::AbstractCell) = isa(cell, CellStart) ? emptyinterval() : (cell.interval, typeof(cell))
isleaf(cell::AbstractCell) = isa(cell, CellEnd) || (isroot(cell) && isempty(cell.children))
AbstractTrees.isroot(cell::AbstractCell) = isa(cell, CellStart)

function BisectableCell(interval::IntervalArithmetic.Interval{T}, parent::Union{CellStart, Cell}) where T<:Number
    return BisectableCell(interval, [], parent)
end

function AtomicCell(interval::IntervalArithmetic.Interval{T}, parent::Union{CellStart, Cell}) where T<:Number
    return AtomicCell(interval, [], parent)
end

function CellStart()
    return CellStart([])
end

function height(cell::AbstractCell)
    if isleaf(cell)
        return 0
    else
        return 1 + height(first(cell.children))
    end
end

function get_constructors(cell::CellStart)
    constructors = []
    current = first(cell.children)
    push!(constructors, constructorof(typeof(current)))
    while !isleaf(current)
        current = first(current.children)
        push!(constructors, constructorof(typeof(current)))
    end
    return constructors
end

import Base: push!

using ConstructionBase

function push!(cell::CellStart, intervals, is_member)
    @assert height(cell) == length(intervals) "Height of the cell must match the number of intervals."
    constructors = get_constructors(cell)
    return push!(cell, intervals, 1, constructors, is_member)
end

function push!(cell::AbstractCell, intervals, i, constructors, is_member)
    l = length(intervals)
    if i > l
        return
    end

    pos = findfirst(isequal(intervals[i]), [child.interval for child in cell.children])
    if !isnothing(pos)
        push!(cell.children[pos], intervals, i+1, constructors, is_member)
    else
        if i == l
            child = constructors[i](intervals[i], cell, is_member)
        else
            child = constructors[i](intervals[i], cell)
            push!(child, intervals, i+1, constructors, is_member)
        end
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
        if isa(cell, BisectableCellEnd)
            return IntervalArithmetic.diam(cell.interval)
        else
            return 0
        end
    elseif isroot(cell)
        return maximum(diam.(cell.children))
    else
        if isa(cell, BisectableCell)
            return maximum([IntervalArithmetic.diam(cell.interval), diam.(cell.children)...])
        else
            return maximum(diam.(cell.children))
        end
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


