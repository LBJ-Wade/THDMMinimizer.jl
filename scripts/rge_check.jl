using THDMMinimizer
using Plots

nvac, cbvac, params = generate_params_and_vacs()
sol = run_parameters(params, params.μ, params.μ + 400.0)

l = @layout [a b c]
p1 = plot(sol.t, sol[1, :], legend=false, framestyle=:box);
p2 = plot(sol.t, sol[2, :], legend=false, framestyle=:box);
p3 = plot(sol.t, sol[3, :], legend=false, framestyle=:box);
plot(p1,p2,p3, layout =l,
     title=[L"m_{11}^2" L"m_{12}^2" L"m_{22}^2"],
     titlefont = font(12))

l = @layout [a b c; d e]
p4 = plot(sol.t, sol[4, :], legend=false, framestyle=:box);
p5 = plot(sol.t, sol[5, :], legend=false, framestyle=:box);
p6 = plot(sol.t, sol[6, :], legend=false, framestyle=:box);
p7 = plot(sol.t, sol[7, :], legend=false, framestyle=:box);
p8 = plot(sol.t, sol[8, :], legend=false, framestyle=:box);
plot(p4,p5,p6,p7,p8, layout =l,
     title=[L"\lambda_{1}" L"\lambda_{2}" L"\lambda_{3}" L"\lambda_{4}" L"\lambda_{5}"],
     titlefont = font(12))

l = @layout [a b c]
p9 = plot(sol.t, sol[9, :], legend=false, framestyle=:box);
p10 = plot(sol.t, sol[10, :], legend=false, framestyle=:box);
p11 = plot(sol.t, sol[11, :], legend=false, framestyle=:box);
plot(p9,p10,p11, layout =l,
     title=[L"y_{t}" L"g'" L"g"],
     titlefont = font(12))
