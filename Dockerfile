FROM julia:latest

WORKDIR /app

# Copy your Julia script into the container
COPY /src /app/src
COPY /examples/ex_5-1_running_example.jl /app/examples/ex_5-1_running_example.jl 
COPY /examples/ex_5-2_stability_controller.jl /app/examples/ex_5-2_stability_controller.jl
COPY /examples/ex_5-3_dubins.jl /app/examples/ex_5-3_dubins.jl
COPY /examples/ex_5-4_circle_collision.jl /app/examples/ex_5-4_circle_collision.jl
COPY /examples/ex_5-5_robot_collision.jl /app/examples/ex_5-5_robot_collision.jl

RUN julia -e 'using Pkg; Pkg.add(name="ArgParse"); Pkg.add(name="IntervalArithmetic", version="0.21.2"); Pkg.add(name="LazySets"); Pkg.add(name="Polyhedra"); Pkg.add(name="StaticArrays"); Pkg.add(name="Symbolics"); Pkg.add(name="CDDLib"); Pkg.add(name="Match"); Pkg.add(name="Plots"); Pkg.add(name="ArgParse"); Pkg.add(name="LaTeXStrings"); Pkg.add(name="BenchmarkTools"); Pkg.add(name="Luxor"); Pkg.add(name="MathTeXEngine");'

# Run the script with arguments
ENTRYPOINT ["julia"]