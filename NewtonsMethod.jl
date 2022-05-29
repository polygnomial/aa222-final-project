"""
Newton's Method

Mean anomaly to eccentric anomaly
"""
function M_to_E(M, e, ε)
    # set up initial eccentric anomaly guess
    if e > 0.5
        E = pi
    else
        E = M
    end
    # set up initial error
    𝛿 = - ((E - e*sin(E) - M) / (1 - e*cos(E)))
    # iterate with Newtons method until reach tolerance
    while(abs(𝛿) > ε)
        𝛿 = - ((E - e*sin(E) - M) / (1 - e*cos(E)))
        E += 𝛿
    end
    return E
end

"""
Calculate true anomaly from time array and starting time
Time array in MDJ  !!!!!
"""
function propogateTrueAnomaly(time_array, orbit::KeplarianOrbit, error = 10^-12)
    # turn time arrat into difference time array in seconds
    Δt = (time_array .- time_array[1]) * 24. * 60. * 60.
    # store orbital values
    e = orbit.e
    a = orbit.a
    ν = orbit.ν
    n = (μ/a^3)^0.5 # mean motion
    # calculate initial mean anomaly
    divider = 1 + e*cos(ν)
    y = (sqrt(1-e^2)*sin(ν))/divider
    x = (e+cos(ν))/divider
    M_0 = atan(y,x) - e*y
    # propogate M for given time array
    M_array = M_0 .+ (Δt .* n)
    # convert time propogated mean anomalies into true anomaly array
    # uses newtons method ...
    nu_array = zeros(0)
    for mean_anomaly in M_array
        E = M_to_E(mod2pi(mean_anomaly), e, error)
        nu = atan((1 - e^2)^0.5 * sin(E), cos(E) - e)
        append!(nu_array, nu)
    end
    return(nu_array)
end
