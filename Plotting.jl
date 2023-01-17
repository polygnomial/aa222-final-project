function convert_orbit_hist_to_ECI(orbit_hist)
    x_hist = []
    y_hist = []
    z_hist = []
    for orbit in orbit_hist
        orbit = EquinoctialOrbit(orbit[1],orbit[2],orbit[3],orbit[4],orbit[5],orbit[6])
        orbit = Equinoctial2ECI(orbit::EquinoctialOrbit)
        push!(x_hist, orbit[1])
        push!(y_hist, orbit[2])
        push!(z_hist, orbit[3])
    end
    return x_hist, y_hist, z_hist
end

function plot_Q_contour_a_e(orbit_initial::KeplarianOrbit, orbit_target::KeplarianOrbit, Q_Params, Sat_Params, orbit_hist)
    a_range = collect(range(orbit_initial.a - 100 * 10^3, orbit_target.a + 100 * 10^3, length = 10^2))
    e_range = collect(range(orbit_initial.e - 0.00001, orbit_target.e + 0.00001, length = 10^2))
    cos_term = cos(orbit_initial.ω + orbit_initial.Ω)
    sin_term = sin(orbit_initial.ω + orbit_initial.Ω)
    target_orbit = Keplarian2Equinoctial(orbit_target)
    target_orbit = [target_orbit.p,
                    target_orbit.f,
                    target_orbit.g,
                    target_orbit.h,
                    target_orbit.k,
                    target_orbit.L]
    orbit = Keplarian2Equinoctial(orbit_initial)
    orbit = [orbit.p,
            orbit.f,
            orbit.g,
            orbit.h,
            orbit.k,
            orbit.L]

    Q_values = []
    for a in a_range
        Q_row = []
        for e in e_range
            # convert keplarian ranges to equinoctial
            p = a * (1 - e^2)
            f = e * cos_term
            g = e * sin_term
            Q_val = log(Q([p, f, g, orbit[4], orbit[5], orbit[6]], target_orbit, Q_Params, Sat_Params))
            push!(Q_row, Q_val)
        end
        push!(Q_values, Q_row)
    end

    a_hist = []
    e_hist = []
    for orbit in orbit_hist
        a = orbit[1]/(1-orbit[2]^2 - orbit[3]^2)
        e = sqrt(orbit[2]^2 + orbit[3]^2)
        push!(a_hist, a)
        push!(e_hist, e)
    end

    return(a_range, e_range, Q_values, a_hist, e_hist)
end
