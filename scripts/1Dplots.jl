using THDMMinimizer
using DelimitedFiles
using LaTeXStrings
import PyPlot; const plt = PyPlot;

data1 = readdlm(string(@__DIR__) * "/../data/verified_a1.csv", ',', skipstart=1);
data2 = readdlm(string(@__DIR__) * "/../data/verified_a2.csv", ',', skipstart=1);

function get_vtree_veffs(ts, datanum, type)
    data = type == :normal ? data2[datanum, :] : data1[datanum, :]

    pars = Params([data[1:12]; 0.0]);
    nvac = Vacuum(data[13], data[14], data[15], 0.0, NotSet);
    cbvac = Vacuum(data[16], data[17], data[18], 0.0, NotSet);
    pars.yt = data[10]

    v1(t) = (1-t) * nvac.v1 + t * cbvac.v1;
    v2(t) = (1-t) * nvac.v2 + t * cbvac.v2;
    α(t) =  (1-t) * nvac.α +  t * cbvac.α;
    fs(t) = Fields([v1(t), v2(t), α(t), 0.0, 0.0, 0.0, 0.0, 0.0])

    vtrees = [potential_tree(fs(t), pars) for t in ts]
    veffs = [potential_eff(fs(t), pars) for t in ts]
    return vtrees, veffs
end

# 16, 22, 26, 31,34, 51, 60, 71, 72, 83, 86, 91, 927, 924, 919, 917, 909, 900, 897
point = 1890
plt.close_figs()
begin
    ts = range(-0.2, stop=1.1, length=100)
    vtrees, veffs = get_vtree_veffs(ts, point, :cb)
    plt.figure(dpi=100)
    plt.xlabel(L"t", fontsize=16)
    plt.ylabel(L"$V(\phi(t))$", fontsize=16)
    plt.plot(ts, vtrees, "--", label=L"$V_{\mathrm{tree}}$", lw=2)
    plt.plot(ts, veffs, label=L"$V_{\mathrm{eff}}$", lw=2)
    plt.xlim([minimum(ts), maximum(ts)])
    plt.grid(true)
    plt.legend(fontsize=14)
    plt.gcf()
    #plt.savefig("cb_1D_" * string(point) * ".pdf")
    #title!(L"$\phi(t) = (1-t)\phi_{\mathrm{EW}} + t\phi_{\mathrm{CB}}$")
    #annotate!([(0, -4.66e8, text(L"EW")), (1.0, -4.75e8, text(L"CB"))])
end

nvac,_,pars=row_to_vacs_and_params(data2[24, :])
sqrt(nvac.v1^2 + nvac.v2^2)
# 2, 9, 16, 19, 24, 26
point = 1042
plt.close_figs()
begin
    ts = range(-0.2, stop=1.1, length=100)
    vtrees, veffs = get_vtree_veffs(ts, point, :normal)
    plt.figure(dpi=100)
    plt.xlabel(L"t", fontsize=16)
    plt.ylabel(L"$V(\phi(t))$", fontsize=16)
    plt.plot(ts, vtrees, "--", label=L"$V_{\mathrm{tree}}$", lw=2)
    plt.plot(ts, veffs, label=L"$V_{\mathrm{eff}}$", lw=2)
    plt.xlim([minimum(ts), maximum(ts)])
    plt.grid(true)
    plt.legend(fontsize=14)
    plt.gcf()
    plt.savefig("normal_1D_" * string(point) * ".pdf")
    #title!(L"$\phi(t) = (1-t)\phi_{\mathrm{EW}} + t\phi_{\mathrm{CB}}$")
    #annotate!([(0, -6.43e8, text(L"EW")), (1.0, -6.38e8, text(L"CB"))])
end


macro horner(x, p...)
    ex = esc(p[end])
    for i = length(p)-1:-1:1
        ex = :(muladd(t, $ex, $(esc(p[i]))))
    end
    Expr(:block, :(t = $(esc(x))), ex)
end

using BenchmarkTools
f(x) = @horner(x, 1, 2, 3, 0, 5)
g(x) = 1 + 2x + 3x^2 + 0*x^3 + 5x^4

@benchmark f(3.5)
@benchmark g(3.5)
