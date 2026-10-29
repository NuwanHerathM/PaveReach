FROM julia:1.12.7

WORKDIR /app

# Copy Julia files into the container
COPY ./src ./src
COPY ./examples/exp_5-1_classification.jl ./examples/exp_5-1_classification.jl 
COPY ./examples/exp_5-1_global_robustness_landscape_class_1.jl ./examples/exp_5-1_global_robustness_landscape_class_1.jl 
COPY ./examples/exp_5-1_global_robustness_landscape_class_2.jl ./examples/exp_5-1_global_robustness_landscape_class_2.jl 
COPY ./examples/exp_5-1_global_robustness_landscape_complete.jl ./examples/exp_5-1_global_robustness_landscape_complete.jl 
COPY ./examples/exp_5-2_classification.jl ./examples/exp_5-2_classification.jl 
COPY ./examples/exp_5-2_global_robustness_landscape.jl ./examples/exp_5-2_global_robustness_landscape.jl 
COPY ./examples/exp_5-3_acasxu_phi_1_1D.jl ./examples/exp_5-3_acasxu_phi_1_1D.jl 
COPY ./examples/exp_5-3_acasxu_phi_1_2D.jl ./examples/exp_5-3_acasxu_phi_1_2D.jl 
COPY ./examples/exp_5-3_acasxu_phi_1_3D.jl ./examples/exp_5-3_acasxu_phi_1_3D.jl 
COPY ./examples/exp_5-3_acasxu_phi_1_4D.jl ./examples/exp_5-3_acasxu_phi_1_4D.jl 
COPY ./examples/exp_5-3_acasxu_clear_distance.jl ./examples/exp_5-3_acasxu_clear_distance.jl 
COPY ./examples/exp_5-3_acasxu_safely_clear_distance.jl ./examples/exp_5-3_acasxu_safely_clear_distance.jl 
COPY ./examples/exp_5-3_acasxu_clear_distance_coarse.jl ./examples/exp_5-3_acasxu_clear_distance_coarse.jl 
COPY ./examples/exp_5-3_acasxu_safely_clear_distance_faster.jl ./examples/exp_5-3_acasxu_safely_clear_distance_faster.jl
# Copy models into the container
COPY ./examples/models/running_example.nnet ./examples/models/running_example.nnet
COPY ./examples/models/normed_etcs_2x250_nodes.onnx ./examples/models/normed_etcs_2x250_nodes.onnx
COPY ./examples/models/ACASXU_run2a_1_1_batch_2000.nnet ./examples/models/ACASXU_run2a_1_1_batch_2000.nnet

# Add Julia packages
RUN julia -e 'using Pkg; Pkg.add(name="ArgParse", version="1.2.0"); Pkg.add(name="IntervalArithmetic", version="0.21.2"); Pkg.add(name="LazySets", version="2.14.2"); Pkg.add(name="Symbolics", version="6.31.0"); Pkg.add(name="ArgParse", version="1.2.0"); Pkg.add(name="BenchmarkTools", version="1.8.0"); Pkg.add(name="Luxor", version="4.3.0"); Pkg.add(name="Match", version="2.4.1"); Pkg.add(name="ONNX", version="0.3.0"); Pkg.add(name="Plots", version="1.41.5"); Pkg.add(name="ReusePatterns", version="0.3.1"); Pkg.add(name="StaticArrays", version="1.9.24"); Pkg.add(name="SoleLogics", version="0.13.7"); Pkg.add(name="ProtoBuf", version="1.2.3"); Pkg.add(url="https://github.com/sisl/NeuralVerification.jl", rev="master");'

# Run the script with arguments
ENTRYPOINT ["julia"]