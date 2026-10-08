function volume_box(box)
    return prod(IntervalArithmetic.diam.(box))
end

function volume_boxes(boxes)
    if isempty(boxes)
        return 0
    end
    return sum(volume_box.(boxes))
end

function merge_intervals(intervals)
    disconnected = deepcopy(intervals)
    merged = []
    while !isempty(disconnected)
        current = pop!(disconnected)
        i = 1
        while i <= length(disconnected)
            if current.lo == disconnected[i].hi
                current = interval(disconnected[i].lo, current.hi)
                popat!(disconnected, i)
                push!(disconnected, current)
                break
            elseif current.hi == disconnected[i].lo
                current = interval(current.lo, disconnected[i].hi)
                popat!(disconnected, i)
                push!(disconnected, current)
                break
            else
                i += 1
            end
        end
        if i > length(disconnected)
            push!(merged, current)
        end
    end
    return merged
end

# Plots drawing

rectangle(p, q) = Shape([p[1],q[1],q[1],p[1]], [p[2],p[2],q[2],q[2]])

function draw_lines(pl, intervals, color)
    for interval in intervals
        p = [interval.lo, -0.1]
        q = [interval.hi, 0.1]
        plot!(pl, rectangle(p,q), color=color, linecolor=nothing, legend=:false)
    end
end

function draw_rows(p, boxes, color)
    draw_lines(p, merge_intervals([box[1] for box in boxes]), color)
end

draw_inn_lines(p, inns) = draw_rows(p, inns, :green)
draw_out_lines(p, outs) = draw_rows(p, outs, :cyan)
draw_delta_lines(p, deltas) = draw_rows(p, deltas, :yellow)

function draw_rectangles(boxes, color)
    for box in boxes
        x_interval = box[1]
        y_interval = box[2]
        p = (x_interval.lo, y_interval.lo)
        q = (x_interval.hi, y_interval.hi)
        plot!(rectangle(p, q), color=color, legend=:false)
    end
end

draw_delta_rectangles(delta) = draw_rectangles(delta, :yellow)
draw_inn_rectangles(inn) = draw_rectangles(inn, :green)
draw_out_rectangles(out) = draw_rectangles(out, :cyan)

function draw(p, X_0, inn, out, delta)
    if isa(X_0, IntervalBox{1, <:Real})
        xticks!((X_0[1].lo:1:X_0[1].hi))
        yaxis!(false)

        draw_delta_lines(p, delta)
        draw_inn_lines(p, inn)
        draw_out_lines(p, out)
    elseif isa(X_0, IntervalBox{2, <:Real})
        xs = X_0[1]
        ys = X_0[2]
        xlims!((xs.lo, xs.hi))
        ylims!((ys.lo, ys.hi))

        draw_delta_rectangles(delta)
        draw_inn_rectangles(inn)
        draw_out_rectangles(out)
    else
        error("Plotting is only supported for 1D and 2D problems.")
    end
end

function print_inn_out_delta(inn, out, delta)
    Base.println("Union of elements of inn: ", merge_intervals([box[1] for box in inn]))
    Base.println("Union of elements of out: ", merge_intervals([box[1] for box in out]))
    Base.println("Union of elements of delta: ", merge_intervals([box[1] for box in delta]))
end

function print_delta_width(delta)
    Base.println("Width of delta regions: ", IntervalArithmetic.diam.(merge_intervals([box[1] for box in delta])))
end

# Luxor drawing

function luxor_box2pq(box)
    x = box[1]
    y = box[2]
    p = Luxor.Point(x.lo, y.lo)
    q = Luxor.Point(x.hi, y.hi)
    return p, q
end

function luxor_rescale(p, q, X_0, width, height, buffer)
    if isa(X_0, IntervalArithmetic.Interval{<:Real})
        scale_x = width / (X_0.hi - X_0.lo)
        scale_y = height / 0.2
        p_rescaled = Luxor.Point(buffer + (p.x - X_0[1].lo) * scale_x, buffer + height - (p.y + 0.1) * scale_y)
        q_rescaled = Luxor.Point(buffer + (q.x - X_0[1].lo) * scale_x, buffer + height - (q.y + 0.1) * scale_y)
    elseif isa(X_0, IntervalBox{1, <:Real})
        scale_x = width / (X_0[1].hi - X_0[1].lo)
        scale_y = height / 0.2
        p_rescaled = Luxor.Point(buffer + (p.x - X_0[1].lo) * scale_x, buffer + height - (p.y + 0.1) * scale_y)
        q_rescaled = Luxor.Point(buffer + (q.x - X_0[1].lo) * scale_x, buffer + height - (q.y + 0.1) * scale_y)
    elseif isa(X_0, IntervalBox{2, <:Real})
        scale_x = width / (X_0[1].hi - X_0[1].lo)
        scale_y = height / (X_0[2].hi - X_0[2].lo)
        p_rescaled = Luxor.Point(buffer + (p.x - X_0[1].lo) * scale_x, buffer + height - (p.y - X_0[2].lo) * scale_y)
        q_rescaled = Luxor.Point(buffer + (q.x - X_0[1].lo) * scale_x, buffer + height - (q.y - X_0[2].lo) * scale_y)
    end
    return p_rescaled, q_rescaled
end

function luxor_rescaled_pq(box, X_0, width, height, buffer)
    p, q = luxor_box2pq(box)
    return luxor_rescale(p, q, X_0, width, height, buffer)
end

function luxor_draw_rows(boxes, color, X_0, width, height, buffer)
    sethue(color)
    for box in boxes
        p = Luxor.Point(box[1].lo, -0.1)
        q = Luxor.Point(box[1].hi, 0.1)
        p_rescaled, q_rescaled = luxor_rescale(p, q, X_0, width, height, buffer)
        Luxor.box(p_rescaled, q_rescaled, :fill)
    end
end

luxor_draw_inn_rows(inn, X_0, width, height, buffer) = luxor_draw_rows(inn, "green", X_0, width, height, buffer)
luxor_draw_out_rows(out, X_0, width, height, buffer) = luxor_draw_rows(out, "cyan", X_0, width, height, buffer)
luxor_draw_delta_rows(delta, X_0, width, height, buffer) = luxor_draw_rows(delta, "yellow", X_0, width, height, buffer)

function luxor_draw_boxes(boxes, color, X_0, width, height, buffer)
    sethue(color)
    for box in boxes
        p, q = luxor_rescaled_pq(box, X_0, width, height, buffer)
        Luxor.box(p, q, :fill)
    end
end

luxor_draw_inn_boxes(inn, X_0, width, height, buffer) = luxor_draw_boxes(inn, "green", X_0, width, height, buffer)
luxor_draw_out_boxes(out, X_0, width, height, buffer) = luxor_draw_boxes(out, "cyan", X_0, width, height, buffer)
luxor_draw_delta_boxes(delta, X_0, width, height, buffer) = luxor_draw_boxes(delta, "yellow", X_0, width, height, buffer)

function luxor_draw(X_0, inn, out, delta, width, height, buffer, xticks=nothing, yticks=nothing)
    if isa(X_0, IntervalArithmetic.Interval{<:Real}) || isa(X_0, IntervalBox{1, <:Real})
        background("white")

        luxor_draw_inn_rows(inn, X_0, width, height, buffer)
        luxor_draw_out_rows(out, X_0, width, height, buffer)
        luxor_draw_delta_rows(delta, X_0, width, height, buffer)

        sethue("black")
        # xticks
        x_start = isnothing(xticks) ? X_0[1].lo : xticks[1]
        x_finish = isnothing(xticks) ? X_0[1].hi : xticks[2]
        tickline(Luxor.Point(buffer, buffer + height), Luxor.Point(buffer + width, buffer + height), startnumber= x_start, finishnumber=x_finish, major=4, minor=0)
    elseif isa(X_0, IntervalBox{2, <:Real})
        background("white")
        
        luxor_draw_inn_boxes(inn, X_0, width, height, buffer)
        luxor_draw_out_boxes(out, X_0, width, height, buffer)
        luxor_draw_delta_boxes(delta, X_0, width, height, buffer)
        
        sethue("black")
        # xticks
        x_start = isnothing(xticks) ? X_0[1].lo : xticks[1]
        x_finish = isnothing(xticks) ? X_0[1].hi : xticks[2]
        tickline(Luxor.Point(buffer, buffer + height), Luxor.Point(buffer + width, buffer + height), startnumber= x_start, finishnumber=x_finish, major=4, minor=0)
        # yticks
        y_start = isnothing(yticks) ? X_0[2].lo : yticks[1]
        y_finish = isnothing(yticks) ? X_0[2].hi : yticks[2]
        tickline(Luxor.Point(buffer + width, buffer + height), Luxor.Point(buffer + width, buffer), startnumber= y_start, finishnumber=y_finish, major=4, minor=0)
    else
        error("Plotting is only supported for 1D and 2D problems.")
    end
end

function save_drawing(X_0::IntervalBox{N, <:Real}, inn, out, delta, filename, xticks=nothing, yticks=nothing) where N
    outdir = "outputs/"
    if !isdir(outdir)
        mkdir(outdir)
        println("Created $(outdir).")
    end
    filepath = outdir * filename
    width = 1000
    height = @match N begin
        1 => 50
        2 => 1000
        _ => error("Plotting is only supported for 1D and 2D problems.")
    end
    buffer = 50
    Drawing(width + 2*buffer, height + 2*buffer, filepath)
    luxor_draw(X_0, inn, out, delta, width, height, buffer, xticks, yticks)
    finish()
    println("The result was saved in $(filepath).")
end

save_drawing(X_0::AbstractVector{IntervalArithmetic.Interval{T}}, inn, out, delta, filename, xticks=nothing, yticks=nothing) where T <: Real = save_drawing(IntervalBox(X_0), inn, out, delta, filename, xticks, yticks)
save_drawing(X_0::IntervalArithmetic.Interval{T}, inn, out, delta, filename, xticks=nothing, yticks=nothing) where T <: Real = save_drawing(IntervalBox(X_0), inn, out, delta, filename, xticks, yticks)

# Logging

function save_log(time, filename)
    outdir = "outputs/"
    if !isdir(outdir)
        mkdir(outdir)
        println("Created $(outdir).")
    end
    filepath = outdir * filename
    open(filepath, "w") do io
        println(io, time)
    end
    println("The time was saved in $(filepath).")
end
