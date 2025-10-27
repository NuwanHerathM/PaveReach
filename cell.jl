using AbstractTrees
using Match

abstract type AbstractCell end

struct OrCell <: AbstractCell
    interval::IntervalArithmetic.Interval{T} where T<:Number
    children::Vector{AbstractCell}
    parent::AbstractCell
end

struct AndCell <: AbstractCell
    interval::IntervalArithmetic.Interval{T} where T<:Number
    children::Vector{AbstractCell}
    parent::AbstractCell
end

mutable struct AtomicCell <: AbstractCell
    const interval::IntervalArithmetic.Interval{T} where T<:Number
    child::Union{AbstractCell, Nothing}
    const parent::AbstractCell
end

struct OrCellStart <: AbstractCell
    children::Vector{AbstractCell}
end

struct AndCellStart <: AbstractCell
    children::Vector{AbstractCell}
end

struct CellEnd <: AbstractCell
    interval::IntervalArithmetic.Interval{T} where T<:Number
    parent::AbstractCell
    is_member::Function
end

struct ConjunctionCell <: AbstractCell
    children::Vector{AbstractCell}
end

struct DisjunctionCell <: AbstractCell
    children::Vector{AbstractCell}
end

ConnectiveCell = Union{ConjunctionCell, DisjunctionCell}

CellStart = Union{OrCellStart, AndCellStart, ConnectiveCell}

AbstractTrees.children(cell::AbstractCell) = isa(cell, CellEnd) ? [] : isa(cell, AtomicCell) ? [cell.child] : cell.children
AbstractTrees.parent(cell::AbstractCell) = isa(cell, CellStart) ? nothing : cell.parent
AbstractTrees.nodevalue(cell::AbstractCell) = @match typeof(cell) begin
    $ConjunctionCell => "∧"
    $DisjunctionCell => "∨"
    $AndCellStart    => (emptyinterval(), typeof(cell))
    $OrCellStart     => (emptyinterval(), typeof(cell))
    _                => (cell.interval, typeof(cell))
    end
isleaf(cell::AbstractCell) = isa(cell, CellEnd)
AbstractTrees.isroot(cell::AbstractCell) = isa(cell, CellStart)

function OrCell(interval::IntervalArithmetic.Interval{T}, parent::AbstractCell) where T<:Number
    return OrCell(interval, [], parent)
end

function AndCell(interval::IntervalArithmetic.Interval{T}, parent::AbstractCell) where T<:Number
    return AndCell(interval, [], parent)
end

function AtomicCell(interval::IntervalArithmetic.Interval{T}, parent::AbstractCell) where T<:Number
    return AtomicCell(interval, nothing, parent)
end

function OrCellStart()
    return OrCellStart([])
end

function AndCellStart()
    return AndCellStart([])
end

function height(cell::AbstractCell)
    if isleaf(cell)
        return 0
    else
        return 1 + height(first(children(cell)))
    end
end

function depth(cell::AbstractCell)
    if isroot(cell)
        return 0
    else
        return 1 + depth(cell.parent)
    end
end

function get_constructors(cell::CellStart)
    constructors = []
    current = first(children(cell))
    push!(constructors, constructorof(typeof(current)))
    while !isleaf(current)
        current = first(children(current))
        push!(constructors, constructorof(typeof(current)))
    end
    return constructors
end

function add_child(cell::AbstractCell, child::AbstractCell)
    push!(cell.children, child)
end

function add_child(cell::AtomicCell, child::AbstractCell)
    @assert isnothing(cell.child) "AtomicCell can only have one child."
    cell.child = child
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

    pos = findfirst(isequal(intervals[i]), [child.interval for child in children(cell)])
    if !isnothing(pos)
        push!(children(cell)[pos], intervals, i+1, constructors, is_member)
    else
        while i < l
            child = constructors[i](intervals[i], cell)
            add_child(cell, child)
            cell = child
            i += 1
        end
        child = constructors[end](intervals[end], cell, is_member)
        add_child(cell, child)
    end
end

function remove!(cell::AbstractCell)
    if isroot(cell)
        return
    end

    siblings = children(cell.parent)
    filter!(sibling -> sibling !== cell, siblings)

    if isempty(siblings)
        remove!(cell.parent)
    end
end

function remove!(cell::CellEnd, intervals)
    if isempty(intervals)
        remove!(cell)
    end
end

function remove!(cell::AbstractCell, intervals)
    if isempty(intervals)
        return
    end
    
    pos = findfirst(isequal(intervals[begin]), [child.interval for child in children(cell)])
    if !isnothing(pos)
        remove!(cell.children[pos], intervals[2:end])
    end
end

function diam(cell::AbstractCell)
    @match typeof(cell) begin
        $CellEnd         => return IntervalArithmetic.diam(cell.interval)
        $AtomicCell      => return maximum(diam.(children(cell)))
        $OrCellStart     => return maximum(diam.(children(cell)))
        $AndCellStart    => return maximum(diam.(children(cell)))
        $ConjunctionCell => return maximum(diam.(children(cell)))
        $DisjunctionCell => return maximum(diam.(children(cell)))
        _             => return maximum([IntervalArithmetic.diam(cell.interval), diam.(children(cell))...])
    end
end

function bisectable_cells(cell::AbstractCell)
    to_visit::Vector{AbstractCell} = [cell]
    cells = []
    while !isempty(to_visit)
        current = pop!(to_visit)
        if isleaf(current)
            continue
        end
        if isa(current, AtomicCell)
            push!(to_visit, current.child)
        else
            append!(to_visit, children(current))
            append!(cells, children(current))
        end
    end
    return cells
end

function largest_bisectable_cell(root::CellStart)
    candidates = bisectable_cells(root)
    diams = IntervalArithmetic.diam.([cell.interval for cell in candidates])
    _, pos = findmax(diams)
    return candidates[pos]
end

function starting_intervals(cell::AbstractCell)
    intervals = []
    while !isa(cell, CellEnd)
        cell_children = children(cell)
        push!(intervals, reduce(∪, [child.interval for child in cell_children]))
        cell = first(cell_children)
    end
    return intervals
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

function make_approximationcell_leaf(box::IntervalBox, inner, outer, parent::AbstractCell)
    cell = ApproximationCell(box, [], parent, inner, outer)
    push!(parent.children, cell)
    return cell
end

#-----------------------------------------------------------------------------------------
