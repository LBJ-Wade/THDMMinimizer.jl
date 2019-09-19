using THDMMinimizer
using DelimitedFiles
import PyPlot; const plt=PyPlot;
using LaTeXStrings
using HomotopyContinuation

data1 = readdlm(string(@__DIR__) * "/../data/verified_a1.csv", ',', skipstart=1);
data2 = readdlm(string(@__DIR__) * "/../data/verified_a2.csv", ',', skipstart=1);

function get_tree_roots(params::Params)
    m112, m122, m222 = params.m112, params.m122, params.m222
    λ1, λ2, λ3, λ4, λ5 = params.λ1, params.λ2, params.λ3, params.λ4, params.λ5

    @polyvar v1 v2 v3

    result = HomotopyContinuation.solve(
        [(1.0 * m112 * v1 -
          1.0 * m122 * v2 +
          0.5 * λ1 * v1 * v1 +
          0.5 * λ1 * v1 * v3 * v3 +
          0.5 * λ3 * v1 * v2 * v2 +
          0.5 * λ4 * v1 * v2 * v2 +
          0.5 * λ5 * v1 * v2 * v2),
         (-1.0 * m122 * v1 +
          1.0 * m222 * v2 +
          0.5 * λ2 * v2 * v2 +
          0.5 * λ3 * v1 * v1 * v2 +
          0.5 * λ3 * v2 * v3 * v3 +
          0.5 * λ4 * v1 * v1 * v2 +
          0.5 * λ5 * v1 * v1 * v2),
         (1.0 * m112 * v3 +
          0.5 * λ1 * v3 * v3 +
          0.5 * λ1 * v1 * v1 * v3 +
          0.5 * λ3 * v2 * v2 * v3)])

    real_solutions(result)
end

function row_to_vacs_and_params(row; tol=1e-5)
    if length(row) == 18
        pars = Params([row[i] for i in 1:12])
        nvac = Vacuum(row[13], row[14], row[15], 0.0, NotSet)
        cbvac = Vacuum(row[16], row[17], row[18], 0.0, NotSet)
        categorize_vacuum!(nvac, pars; tol=tol)
        categorize_vacuum!(cbvac, pars; tol=tol)
        nvac.potential = potential_eff(Fields(nvac), pars)
        cbvac.potential = potential_eff(Fields(cbvac), pars)
        return nvac, cbvac, pars
    elseif length(row) == 15
        pars = Params([row[i] for i in 1:12])
        vac = Vacuum(row[13], row[14], row[15], 0.0, NotSet)
        categorize_vacuum!(vac, pars; tol=tol)
        vac.potential = potential_eff(Fields(vac), pars)
        return vac, pars
    else
        throw("row doesn't contain correct number of columns")
    end
end

function get_vtree_veffs_plane(ts, ss, θ, data)

    nvac, cbvac, pars = row_to_vacs_and_params(data)
    #tree_roots = get_tree_roots(params)

    X, Y, Z = -(nvac.v1 - cbvac.v1), -(nvac.v2 - cbvac.v2), -(nvac.α - cbvac.α)
    mag = sqrt(X^2 + Y^2 + Z^2)
    x = X / mag
    y = Y / mag
    z = Z / mag
    sx = y*cos(θ) + sqrt(1 - x^2 - y^2)*sin(θ)
    sy = ((-1 - x + y^2)*cos(θ) + y*sqrt(1 - x^2 - y^2)*sin(θ))/(1 + x)
    sz = (y*sqrt(1 - x^2 - y^2)*cos(θ) - (x + x^2 + y^2)*sin(θ))/(1 + x)

    v1(t, s) = (t * x + sx * s) * mag + nvac.v1;
    v2(t, s) = (t * y + sy * s) * mag + nvac.v2;
    α(t, s) =  (t * z + sz * s) * mag + nvac.α;
    fs(t, s) = Fields([v1(t, s), v2(t, s), α(t, s), 0.0, 0.0, 0.0, 0.0, 0.0])
    vtrees = [potential_tree(fs(t, s), pars) for t in ts, s in ss]
    veffs = [potential_eff(fs(t, s), pars) for t in ts, s in ss]
    tree_roots =
    return vtrees, veffs, potential_eff(fs(0.0, 0.), pars), potential_eff(fs(1.0, 0.0), pars)
end

point=897
plt.close_figs()
begin
    plt.figure(dpi=100, figsize=(9,6))
    plt.suptitle(L"$V_{\mathrm{eff}}(\phi_{CB})<V_{\mathrm{eff}}(\phi_{EW})$", fontsize=16)
    plt.subplot(121)
    plt.title("Effective Potential", fontsize=16)
    ts = range(-0.2, stop=1.2, length=100)
    ss = range(-0.2, stop=0.2, length=100)
    vtrees, veffs, vn, vcb = get_vtree_veffs_plane(ts, ss, 0, data1[point,:])
    plt.contourf(ts, ss, veffs', cmap="viridis", levels=10)
    plt.contour(ts, ss, veffs', cmap="binary",levels=35, linewidths=1,alpha=0.4)
    plt.scatter([0.0], [0.0], c="white")
    plt.scatter([1.0], [0.0], c="white")
    plt.text(0.90, 0.02, L"\mathrm{CB}", fontsize=30, fontdict=Dict(:color=>"white"))
    plt.text(-0.1, 0.02, L"\mathrm{EW}", fontsize=30, fontdict=Dict(:color=>"white"))

    plt.subplot(122)
    plt.yticks([])
    plt.title("Tree-Level Potential", fontsize=16)
    plt.contourf(ts, ss, vtrees', cmap="viridis", levels=10, alpha=1.0)
    plt.colorbar()
    plt.contour(ts, ss, vtrees', cmap="Greys", levels=30, linewidths=1,alpha=0.4)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.gcf()
    plt.savefig("cb_2D_" * string(point) * ".pdf")
end




point=26
plt.close_figs()
begin
    plt.figure(dpi=100, figsize=(9,6))
    plt.suptitle(L"$V_{\mathrm{eff}}(\phi_{EW})<V_{\mathrm{eff}}(\phi_{CB})$", fontsize=16)
    plt.subplot(121)
    plt.title("Effective Potential", fontsize=16)
    ts = range(-0.2, stop=1.2, length=100)
    ss = range(-0.4, stop=0.4, length=100)
    vtrees, veffs, vn, vcb = get_vtree_veffs_plane(ts, ss, 0.0, data2[point,:])
    plt.contourf(ts, ss, veffs', cmap="viridis", levels=20)
    plt.contour(ts, ss, veffs', cmap="Greys", levels=60, linewidths=1,alpha=0.4)
    plt.scatter([0.0], [0.0], c="white")
    plt.scatter([1.0], [0.0], c="white")
    plt.text(0.90, 0.02, L"\mathrm{CB}", fontsize=30, fontdict=Dict(:color=>"white"))
    plt.text(-0.1, 0.02, L"\mathrm{EW}", fontsize=30, fontdict=Dict(:color=>"white"))

    plt.subplot(122)
    plt.yticks([])
    plt.title("Tree-Level Potential", fontsize=16)
    plt.contourf(ts, ss, vtrees', cmap="viridis", levels=10)
    plt.colorbar()
    plt.contour(ts, ss, vtrees', cmap="Greys", levels=35, linewidths=1,alpha=0.4)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.gcf()
    plt.savefig("normal_2D_" * string(point) * ".pdf")
end
