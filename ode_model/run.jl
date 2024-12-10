using DifferentialEquations
using Plots

include("parameters.jl")
include("ode_system.jl")
include("protein_spacing.jl")
include("ode_call.jl")

#Number of binding sites / particle
# nps  = collect(range(5, 250, length = 36))
nps = [5, 6, 10, 20, 22, 25, 30, 35,40, 50, 60, 70, 100, 120, 150, 200, 250]
nps = Int.(nps)

values =  [] #Store relative C3a concentration
values2 = [] #Store relative bound C3b concentration
values3 =  [] #Store C3a concentration
values4 = [] #Store bound C3b concentration

tend = 10 #minutes
tspan = (0.0, tend)

for sites_np in nps
    conc_c3b_sites = sites_np*conc_np
    scaling_fact = c_c3b/conc_c3b_sites
    scaling_fact2 = 1.61e-5/conc_c3b_sites
    sol = ode_call(conc_c3b_sites, scaling_fact, scaling_fact2)
    push!(values, (sol[89, end])/sol[10, 1])
    push!(values2, (sol[11, end])/sol[10, 1])    
    push!(values3, (sol[89, end]))
    push!(values4, (sol[11, end]))    
end

Plots.scatter(distance, values, xflip=true, label=false, c = "blue")
display(Plots.plot!(distance, values, xlabel="Protein-Protein spacing (nm)", ylabel="C3a (μM)/Binding Site (μM)",
                    xtickfontsize=14,ytickfontsize=14,yguidefontsize=14,xguidefontsize=14,legendfontsize=14, xflip=true, label=false, c = "blue"))

Plots.scatter(distance, values2, xflip=true, label=false, c = "blue")
display(Plots.plot!(distance, values2, xlabel="Protein-Protein spacing (nm)", ylabel="C3b (μM)/Binding Site (μM)",
                    xtickfontsize=14,ytickfontsize=14,yguidefontsize=14,xguidefontsize=14,legendfontsize=14, xflip=true, label=false, c = "blue"))

Plots.scatter(distance, values3, xflip=true, label=false, c = "blue")
display(Plots.plot!(distance, values3, xlabel="Protein-Protein spacing (nm)", ylabel="C3a (μM)",
                    xtickfontsize=14,ytickfontsize=14,yguidefontsize=14,xguidefontsize=14,legendfontsize=14, xflip=true, label=false, c = "blue"))

Plots.scatter(distance, values4, xflip=true, label=false, c = "blue")
display(Plots.plot!(distance, values4, xlabel="Protein-Protein spacing (nm)", ylabel="C3b (μM)",
                    xtickfontsize=14,ytickfontsize=14,yguidefontsize=14,xguidefontsize=14,legendfontsize=14, xflip=true, label=false, c = "blue"))
