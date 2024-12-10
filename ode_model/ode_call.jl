function ode_call(conc_c3b_sites, scaling_fact, scaling_fact2)
    
    y0 = zeros(107)

    y0[88], y0[98], y0[101], y0[102], y0[103],
    y0[104], y0[90], y0[91], y0[92], y0[93], y0[94],
    y0[95], y0[96], y0[97], y0[105],  y0[106] ,
    y0[107],  y0[10],  y0[50]  = (5.4, 0.37, 0.54, 0.50, 0.36, 0.90, 2.2, 0.083,
                                 0.4, 0.47, 0.0009, 3.2, 0.0083, 0.027, 6.0, 
                                 0.43, 0.21, conc_c3b_sites*1e6, 5.0)    

    prob = ODEProblem(ODEfunc!, y0, tspan, [scaling_fact, scaling_fact2])
    sol = solve(prob, AutoTsit5(Rosenbrock23()))

    return sol
end