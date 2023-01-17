__precompile__()
module JQlaw
    include("Qlaw.jl")
    include("Integrator.jl")
    include("NewtonsMethod.jl")
    include("plotting.jl")

    # Qlaw
    export µ
    export R_e
    export EquinoctialOrbit
    export KeplarianOrbit
    export Equinoctial2Keplarian
    export Equinoctial2ECI
    export Keplarian2Equinoctial
    export GaussVariationalEquationsEquinoctial
    export QParams
    export sat_params
    export propogateTrueAnomaly
    export Qlaw
    export Equinoctial2ECI
    export convert_orbit_hist_to_ECI
    export plot_Q_contour_a_e
    export make_make_objective

    # integrator
    export RK45
end
