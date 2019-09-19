using THDMMinimizer
using DataFrames
using StatsPlots
using CSV
using LaTeXStrings

df_a1 = DataFrame(CSV.read(string(@__DIR__) * "/../data/verified_a1.csv"));
df_a2 = DataFrame(CSV.read(string(@__DIR__) * "/../data/verified_a2.csv"));
df_b = DataFrame(CSV.read(string(@__DIR__) * "/../data/verified_b.csv"));
df_c = DataFrame(CSV.read(string(@__DIR__) * "/../data/verified_c.csv"));

lm112 = L"$m_{11}^2$";
lm122 = L"$m_{12}^2$";
lm222 = L"$m_{22}^2$";
lλ1 = L"$\lambda_{1}$";
lλ2 = L"$\lambda_{2}$";
lλ3 = L"$\lambda_{3}$";
lλ4 = L"$\lambda_{4}$";
lλ5 = L"$\lambda_{5}$";

@df df_a1 corrplot(cols(1:8), grid = false, label=[lm112,lm122,lm222,lλ1,lλ2,lλ3,lλ4,lλ5])
@df df_a2 corrplot(cols(1:8), grid = false, label=[lm112,lm122,lm222,lλ1,lλ2,lλ3,lλ4,lλ5])
@df df_b corrplot(cols(1:8), grid = false, label=[lm112,lm122,lm222,lλ1,lλ2,lλ3,lλ4,lλ5])
@df df_c corrplot(cols(1:8), grid = false, label=[lm112,lm122,lm222,lλ1,lλ2,lλ3,lλ4,lλ5])


row = dfa1[end,:]

pars = Params([row[i] for i in 1:12])
nvac = Vacuum(row[13], row[14], row[15], 0.0, NotSet)
cbvac = Vacuum(row[16], row[17], row[18], 0.0, NotSet)

autodiff_gradient(potential_eff, Fields(nvac), pars, (R1,R2,C1,C2,C3,C4,I1,I2))
autodiff_gradient(potential_eff, Fields(cbvac), pars, (R1,R2,C1,C2,C3,C4,I1,I2))
one_loop_masses(nvac, pars)
one_loop_masses(cbvac, pars)
