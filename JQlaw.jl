__precompile__()
module JQlaw
    include("Qlaw.jl")
    include("Integrator.jl")
    include("NewtonsMethod.jl")

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
    export make_objective

    # integrator
    export RK45
end
