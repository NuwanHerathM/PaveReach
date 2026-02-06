FROM julia:1.12.4

WORKDIR /app

# Copy Julia script into the container
COPY /src /app/src
COPY /examples/ex_5-1_running_example.jl /app/examples/ex_5-1_running_example.jl 
COPY /examples/ex_5-1_running_example_no_timing.jl /app/examples/ex_5-1_running_example_no_timing.jl 
COPY /examples/ex_5-2_stability_controller.jl /app/examples/ex_5-2_stability_controller.jl
COPY /examples/ex_5-3_dubins.jl /app/examples/ex_5-3_dubins.jl
COPY /examples/ex_5-4_circle_collision.jl /app/examples/ex_5-4_circle_collision.jl
COPY /examples/ex_5-5_robot_collision.jl /app/examples/ex_5-5_robot_collision.jl

# Add Julia packages
RUN julia -e 'using Pkg; Pkg.add(name="ArgParse", version="1.2.0"); Pkg.add(name="IntervalArithmetic", version="0.21.2"); Pkg.add(name="LazySets", version="2.14.2"); Pkg.add(name="Symbolics", version="6.31.0"); Pkg.add(name="Match", version="2.4.1"); Pkg.add(name="Plots", version="1.41.5"); Pkg.add(name="LaTeXStrings", version="1.4.0"); Pkg.add(name="BenchmarkTools", version="1.6.3"); Pkg.add(name="Luxor", version="4.3.0"); Pkg.add(name="MathTeXEngine", version="0.6.7");'

# Run the script with arguments
ENTRYPOINT ["julia"]