using ONNX
using NeuralVerification

# Thanks to Francesc for this function!
function read_onnx_mlp(path::AbstractString)
    model = open(path) do io
        ONNX.ProtoBuf.decode(ONNX.ProtoBuf.ProtoDecoder(io), ONNX.ModelProto)
    end

    init = Dict(init.name => ONNX.array(init) for init in model.graph.initializer)
    layers = NeuralVerification.Layer[]
    pending_W = nothing
    pending_b = nothing
    pending_out = nothing

    for node in model.graph.node
        op = String(node.op_type)

        if op == "MatMul"
            pending_W === nothing || error("Unexpected MatMul before finalizing previous layer.")
            weight_names = [name for name in node.input if haskey(init, name)]
            length(weight_names) == 1 || error("MatMul should have exactly one initializer input.")
            W = init[weight_names[1]]
            ndims(W) == 2 || error("MatMul weight must be 2D, got size $(size(W)).")
            pending_W = Matrix{Float64}(W)
            pending_b = nothing
            pending_out = node.output[1]
            continue
        end

        if op == "Gemm"
            pending_W === nothing || error("Unexpected Gemm before finalizing previous layer.")
            attrs = Dict(node.attribute)
            transA = get(attrs, :transA, 0)
            transB = get(attrs, :transB, 0)
            alpha = get(attrs, :alpha, 1.0)
            beta = get(attrs, :beta, 1.0)
            transA == 0 || error("Gemm with transA=$(transA) is not supported.")

            W = init[node.input[2]]
            ndims(W) == 2 || error("Gemm weight must be 2D, got size $(size(W)).")
            W = Matrix{Float64}(W)
            if transB == 1
                W = permutedims(W)
            elseif transB != 0
                error("Gemm with transB=$(transB) is not supported.")
            end
            W .*= alpha
            pending_W = W

            if length(node.input) >= 3 && haskey(init, node.input[3])
                b = Vector{Float64}(vec(init[node.input[3]]))
                pending_b = b .* beta
            else
                pending_b = nothing
            end
            pending_out = node.output[1]
            continue
        end

        if op == "Add" && pending_W !== nothing && pending_out !== nothing && (pending_out in node.input)
            other = node.input[1] == pending_out ? node.input[2] : node.input[1]
            haskey(init, other) || error("Add after MatMul/Gemm must use a bias initializer.")
            pending_b = Vector{Float64}(vec(init[other]))
            pending_out = node.output[1]
            continue
        end

        if op in ("Relu", "Sigmoid", "Tanh", "Identity") && pending_out !== nothing && (pending_out in node.input)
            act = if op == "Relu"
                NeuralVerification.ReLU()
            elseif op == "Sigmoid"
                NeuralVerification.Sigmoid()
            elseif op == "Tanh"
                NeuralVerification.Tanh()
            else
                NeuralVerification.Id()
            end

            pending_W === nothing && error("Activation $(op) without a preceding linear layer.")
            if pending_b === nothing
                pending_b = zeros(Float64, size(pending_W, 1))
            end
            length(pending_b) == size(pending_W, 1) || error("Bias length does not match weight rows.")
            push!(layers, NeuralVerification.Layer(pending_W, pending_b, act))
            pending_W = nothing
            pending_b = nothing
            pending_out = nothing
        end
    end

    if pending_W !== nothing
        if pending_b === nothing
            pending_b = zeros(Float64, size(pending_W, 1))
        end
        length(pending_b) == size(pending_W, 1) || error("Bias length does not match weight rows.")
        push!(layers, NeuralVerification.Layer(pending_W, pending_b, NeuralVerification.Id()))
    end

    return Network(layers)
end
