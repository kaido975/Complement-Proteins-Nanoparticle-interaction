# this script installs all required packages
using Pkg
Pkg.update()
if ! in("DifferentialEquations",keys(Pkg.dependencies())) Pkg.add("DifferentialEquations") end
if ! in("Plots",keys(Pkg.dependencies())) Pkg.add("Plots") end
@warn("Packages installed!")