module THDMMinimizer

using Reexport
using ForwardDiff: Dual, partials, value
using LinearAlgebra
using StatsBase
using NLsolve
using DifferentialEquations
using DelimitedFiles
@reexport using Optim

include("constants.jl")
include("types.jl")
include("jacobi.jl")
include("scalar_squared_masses.jl")
include("fermion_squared_masses.jl")
include("gauge_squared_masses.jl")
include("autodiff.jl")
include("potentials.jl")
include("checks.jl")
include("root_finding.jl")
include("minimizing.jl")
include("init.jl")
include("rge.jl")
include("scanning_tools.jl")

# contatns
export M_TOP
export U1Y_COUP
export SU2_COUP
export HIGGS_VEV
export ZERO_MASS_TOL
export ZERO_DERIV_TOL
# automatic differentiation
export autodiff_gradient
export autodiff_hessian
# Checks
export check_gradients
export check_one_loop_masses
export check_are_perturbative
export categorize_vacuum!
# types
export Fields
export FieldIndex
export R1
export R2
export C1
export C2
export C3
export C4
export I1
export I2
export Params
export ParamIndex
export M112
export M122
export M222
export Λ1
export Λ2
export Λ3
export Λ4
export Λ5
export is_bounded
export Vacuum
export ExtremaType
export Minimum
export Maximum
export Saddle
export Undefined
export NotSet
export set_top_yukawa!
export generate_normal_vac
export generate_cb_vac
# jacobi algorithm
export jacobi
export JacobiMaxIterError
# Potentials
export potential_tree
export potential_eff
export one_loop_masses
# Masses
export scalar_squared_mass_matrix
export scalar_squared_masses
export fermion_squared_masses
export gauge_sqaured_mass_matrix
export gauge_squared_masses
# Solving
export try_find_root_eff!
export minimize!
export solve_tree_eqns!
export generate_params_and_vacs
export solve_root_eqns
export find_new_minima
export catagorize_results
export find_deepest_normal_min
export find_deepest_cb_min
export write_results_to_file
# RGE
export run_parameters



end # module
