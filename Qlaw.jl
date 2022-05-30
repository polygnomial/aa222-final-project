using StaticArrays
using LinearAlgebra

# constants
const µ = 3.986004419 * 10^14 # standard gravitational parameter earth [m^3/s^2]
const R_e = 6371000. # average earth radius [m]
const g_0 = 9.80665  # m/s^2 - standard gravity (for Isp calcs)

struct EquinoctialOrbit
    p # semi latus rectum [meters]
    f
    g
    h
    k
    L # true longitude [rad]
end

struct KeplarianOrbit
    a # semi-major axis [meters]
    e # eccentricity [-]
    i # inclinations [rad]
    ω # arg. of perigee [rad]
    Ω # raan [rad]
    ν # true anomaly [rad]
end

# https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
function Equinoctial2Keplarian(orbit::EquinoctialOrbit)
    # unload orbit parameters
    p = orbit.p
    f = orbit.f
    g = orbit.g
    h = orbit.h
    k = orbit.k
    L = orbit.L
    # convert orbital elements
    a = p/(1-f^2-g^2)
    e = sqrt(f^2 + g^2)
    i = atan(2*sqrt(h^2+k^2), 1-h^2-k^2)
    ω = atan(g*h-f*k, f*h+g*k)
    Ω = atan(k,h)
    ν = L - atan(g,f)
    return(KeplarianOrbit(a, e, i, ω, Ω, ν))
end

function Equinoctial2ECI(orbit::EquinoctialOrbit)
    # unload orbit parameters
    p = orbit.p
    f = orbit.f
    g = orbit.g
    h = orbit.h
    k = orbit.k
    L = orbit.L
    # calculate important values
    α_2 = h^2 - k^2
    s_2 = 1 + h^2 + k^2
    w = 1 + f*cos(L) + g*sin(L)
    r = p/w
    # convert orbital elements
    r_x = (r/s_2)*(cos(L)+α_2*cos(L)+2*h*k*sin(L))
    r_y = (r/s_2)*(sin(L)+α_2*sin(L)+2*h*k*cos(L))
    r_z = (2*r/s_2)*(h*sin(L)-k*cos(L))
    v_x = (-1/s_2)*sqrt(µ/p)*(sin(L)+α_2*sin(L)-2*h*k*cos(L)+g-2*f*h*k+α_2*g)
    v_y = (-1/s_2)*sqrt(µ/p)*(-cos(L)+α_2*cos(L)+2*h*k*sin(L)-f+2*g*h*k+α_2*f)
    v_z = (2/s_2)*sqrt(µ/p)*(h*cos(L)+k*sin(L)+f*h+g*k)
    return([r_x, r_y, r_z, v_x, v_y, v_z])
end

function Keplarian2Equinoctial(orbit::KeplarianOrbit)
    # unload orbit parameters
    a = orbit.a
    e = orbit.e
    i = orbit.i
    ω = orbit.ω
    Ω = orbit.Ω
    ν = orbit.ν
    # convert orbital elements
    p = a*(1-e^2)
    f = e*cos(ω+Ω)
    g = e*sin(ω+Ω)
    h = tan(i/2)*cos(Ω)
    k = tan(i/2)*sin(Ω)
    L = Ω + ω + ν
    return(EquinoctialOrbit(p,f,g,h,k,L))
end

function GaussVariationalEquationsEquinoctial(orbit, time; accel=[0.,0.,0.])
    # unpack orbit
    p = orbit[1] # semi latus rectum [meters]
    f = orbit[2]
    g = orbit[3]
    h = orbit[4]
    k = orbit[5]
    L = orbit[6] # true longitude [rad]

    # calcualte useful quantities
    cosL = cos(L)
    sinL = sin(L)
    w = 1 + f*cosL+g*sinL
    root_p_µ = sqrt(p/µ)
    s_2 = 1 + h^2 + k^2

    # calculate matrix
    A = [0 2*p/w 0;
        sinL (1/w)*((w+1)*cosL+f) -(g/w)*(h*sinL-k*cosL);
        -cosL (w+1)*sinL+g (f/w)*(h*sinL-k*cosL);
        0 0 (s_2*cosL)/(2*w);
        0 0 (s_2*sinL)/(2*w);
        0 0 h*sinL-k*cosL]

    A = root_p_µ * A
    b = [0, 0, 0, 0, 0, sqrt(µ*p)*(w/p)^2]

    # calculate doe/dt
    dOEdt = A*accel + b
end

struct QParams
    W_p    # weight of the minimum radius function
    k      # internal weighting for minimum pariapsis
    rp_min # minimum periapsis radius (penalty term)
    W_a    # Q weight for semi-major axis (0 for free variable)
    W_f    # Q weight for f (0 for free variable)
    W_g    # Q weight for g (0 for free variable)
    W_h    # Q weight for h (0 for free variable)
    W_k    # Q weight for k (0 for free variable)
    m      # scaling weight to keep a from going to inf (nominal value = 3)
    n      # scaling weight to keep a from going to inf (nominal value = 4)
    r      # scaling weight to keep a from going to inf (nominal value = 2)
    central_difference_step # central difference stepsize
    # for following thresholds, 0 would mean always thrust, 1 would mean never thrust
    ƞ_a    # absolute effectivity threshold
    ƞ_r    # relative effectivity threshold
    num_sample_points_per_orbit    # number of sample points to calculate over an orbit to determine effectivity
    step_size # [sec]
    max_iter # max number of integration steps
    Q_convergence # [sec], q convergence threshold
end

struct sat_params
    mass        # mass of satellite [kg]
    prop_ratio  # ratio of prop mass to total mass
    accel       # acceleration [m/s^2]
    isp         # specific impulse [s]
end

function mass_loss(mass, sat_params::sat_params, ∆t)
    f = mass * sat_params.accel
    mdot = f / sat_params.isp / g_0
    return(mass - (mdot * ∆t))
end

function Q(orbit, target_orbit, Q_params, sat_params)
    # pull current orbit and target elements
    p = orbit[1]
    f = orbit[2]
    g = orbit[3]
    h = orbit[4]
    k = orbit[5]
    p_targ = target_orbit[1]
    f_targ = target_orbit[2]
    g_targ = target_orbit[3]
    h_targ = target_orbit[4]
    k_targ = target_orbit[5]
    a = p/(1-f^2-g^2)
    a_targ = p_targ/(1-f_targ^2-g_targ^2)
    # calculate needed intermediate orbital values
    e = sqrt(f^2 + g^2)
    r_p = p/(1+e) # radius pariapsis
    s_2 = 1 + h^2 + k^2
    # calculate Q equation values
    thrust = sat_params.accel
    P = exp(Q_params.k * (1 - (r_p/Q_params.rp_min)))
    S_a = (1 + (abs(a - a_targ)/(Q_params.m * a_targ))^Q_params.n)^(1/Q_params.r)
    a_dot_xx = 2*thrust*a*sqrt((a/µ)*((1+sqrt(f^2 + g^2))/(1-sqrt(f^2 + g^2))))
    f_dot_xx = 2*thrust*sqrt(p/µ)
    g_dot_xx = f_dot_xx
    h_dot_xx = 0.25 * f_dot_xx * (s_2/(sqrt(1-g^2)+f))
    k_dot_xx = 0.25 * f_dot_xx * (s_2/(sqrt(1-f^2)+g))
    # calculate Q
    Q_a = S_a * Q_params.W_a * ((a-a_targ)/a_dot_xx)^2
    Q_f = Q_params.W_f * ((f-f_targ)/f_dot_xx)^2
    Q_g = Q_params.W_g * ((g-g_targ)/g_dot_xx)^2
    Q_h = Q_params.W_h * ((h-h_targ)/h_dot_xx)^2
    Q_k = Q_params.W_k * ((k-k_targ)/k_dot_xx)^2
    Q = (1+(Q_params.W_p*P))*(Q_a+Q_f+Q_g+Q_h+Q_k)
    return Q
end

function calc_D(orbit, target_orbit, Q_params, sat_params)
    # pull current orbit and target elements
    p = orbit[1]
    f = orbit[2]
    g = orbit[3]
    h = orbit[4]
    k = orbit[5]
    L = orbit[6]
    p_targ = target_orbit[1]
    f_targ = target_orbit[2]
    g_targ = target_orbit[3]
    h_targ = target_orbit[4]
    k_targ = target_orbit[5]
    a = p/(1-f^2-g^2)
    a_targ = p_targ/(1-f_targ^2-g_targ^2)
    # calculate needed intermediate orbital values
    e = sqrt(f^2 + g^2)
    n = sqrt(µ/a^3)
    b = a*sqrt(1-e^2)
    h_mom = n*a*b
    w = 1+f*cos(L)+g*sin(L)
    r = p/w
    s_2 = 1 + h^2 + k^2

    # calculate dQ_dae
    Q_of_a(a) = Q([a, orbit[2], orbit[3], orbit[4], orbit[5], orbit[6]],
                  target_orbit,
                  Q_params,
                  sat_params)
    Q_of_f(f) = Q([orbit[1], f, orbit[3], orbit[4], orbit[5], orbit[6]],
                target_orbit,
                Q_params,
                sat_params)
    Q_of_g(g) = Q([orbit[1], orbit[2], g, orbit[4], orbit[5], orbit[6]],
                target_orbit,
                Q_params,
                sat_params)
    Q_of_h(h) = Q([orbit[1], orbit[2], orbit[3], h, orbit[5], orbit[6]],
                target_orbit,
                Q_params,
                sat_params)
    Q_of_k(k) = Q([orbit[1], orbit[2], orbit[3], orbit[4], k, orbit[6]],
                target_orbit,
                Q_params,
                sat_params)
    dQ_da = central_difference(orbit[1], Q_of_a, Q_params.central_difference_step)
    dQ_df = central_difference(orbit[2], Q_of_f, Q_params.central_difference_step)
    dQ_dg = central_difference(orbit[3], Q_of_g, Q_params.central_difference_step)
    dQ_dh = central_difference(orbit[4], Q_of_h, Q_params.central_difference_step)
    dQ_dk = central_difference(orbit[5], Q_of_k, Q_params.central_difference_step)

    dadt_dFr = ((2*a^2)/(h_mom)) * (f*sin(L)-g*cos(L))
    dadt_dFt = ((2*a^2)/(h_mom)) * (p/r)
    dadt_dFn = 0

    dfdt_dFr = sqrt(p/µ)*sin(L)
    dfdt_dFt = sqrt(p/µ)*(1/w)*((w+1)*cos(L)+f)
    dfdt_dFn = -sqrt(p/µ)*(g/w)*(h*sin(L)-k*cos(L))

    dgdt_dFr = -sqrt(p/µ)*cos(L)
    dgdt_dFt = sqrt(p/µ)*(1/w)*((w+1)*sin(L)+g)
    dgdt_dFn = sqrt(p/µ)*(g/w)*(h*sin(L)-k*cos(L))

    dhdt_dFr = 0
    dhdt_dFt = 0
    dhdt_dFn = sqrt(p/µ)*(s_2/(2*w))*cos(L)

    dkdt_dFr = 0
    dkdt_dFt = 0
    dkdt_dFn = sqrt(p/µ)*(s_2/(2*w))*sin(L)

    # calculate D terms
    D1 = dQ_da*dadt_dFt + dQ_df*dfdt_dFt + dQ_dg*dgdt_dFt
    D2 = dQ_da*dadt_dFr + dQ_df*dfdt_dFr + dQ_dg*dgdt_dFr
    D3 = dQ_df*dfdt_dFn + dQ_dg*dgdt_dFn + dQ_dh*dhdt_dFn + dQ_dk*dkdt_dFn

    return(D1, D2, D3)
end

function evaluate_orbit_location(orbit, target_orbit, Q_params, sat_params)
    num_eval_points = Q_params.num_sample_points_per_orbit

    # pull current orbit true longitude
    f = orbit[2]
    g = orbit[3]
    L_current = orbit[6]

    # create array of true longitude to search
    L_array = collect(range(0, (2*pi - (2*pi/num_eval_points)), length = (num_eval_points - 1)))
    push!(L_array, L_current)

    # loop through true longitudes to find dQdt min and max
    D1_array = []
    D2_array = []
    D3_array = []
    for L in L_array
        D1, D2, D3 = calc_D([orbit[1], orbit[2], orbit[3], orbit[4], orbit[5], L],
                            target_orbit, Q_params, sat_params)
        push!(D1_array, D1)
        push!(D2_array, D2)
        push!(D3_array, D3)
    end

    # calculate dQdt min over alpha and beta
    dQdt = []
    for i in 1:length(D1_array)
        dQdt_min = -sqrt(D1_array[i]^2 + D2_array[i]^2 + D3_array[i]^2)
        push!(dQdt, dQdt_min)
    end

    # determine dqdt min and max over alpha beta and true longitude
    dQdt_min = minimum(dQdt)
    dQdt_max = maximum(dQdt)
    dQdt_current = dQdt[end]

    # calculate rlative and absolute efficiencies
    ƞ_a = dQdt_current/dQdt_min
    ƞ_r = (dQdt_current-dQdt_max)/(dQdt_min-dQdt_max)
    if isnan(ƞ_r)
        ƞ_r = 1.
    end

    # check with Q param effetivity thresholds
    if (ƞ_a >= Q_params.ƞ_a) && (ƞ_r >= Q_params.ƞ_r) && (dQdt_current < 0.)
        # thrust in optimal direction for step size
        α_optimal = atan(-D2_array[end], -D1_array[end])
        β_optimal = atan(-D3_array[end]/sqrt(D1_array[end]^2 + D2_array[end]^2))
        F_t = sat_params.accel * cos(β_optimal) * cos(α_optimal)
        F_r = sat_params.accel * cos(β_optimal) * sin(α_optimal)
        F_n = sat_params.accel * sin(β_optimal)
        return (true, F_t, F_r, F_n)
    else
        # do not thrust, dwell for step size
        return (false, 0, 0, 0)
    end
end

function central_difference(input, function_name, step_size)
    derivative = (function_name(input+step_size) - function_name(input-step_size))/(2*step_size)
    return(derivative)
end

function Qlaw(orbit_inital::KeplarianOrbit, orbit_target::KeplarianOrbit, Q_Params::QParams, Sat_Params::sat_params, time_initial)
    # convert keplarian orbits into equinoctial orbits, store orbits as arrays
    orbit_initial_equinoctial = Keplarian2Equinoctial(orbit_inital)
    orbit_target_equinoctial = Keplarian2Equinoctial(orbit_target)

    orbit = [orbit_initial_equinoctial.p,
            orbit_initial_equinoctial.f,
            orbit_initial_equinoctial.g,
            orbit_initial_equinoctial.h,
            orbit_initial_equinoctial.k,
            orbit_initial_equinoctial.L]

    orbit_target = [orbit_target_equinoctial.p,
            orbit_target_equinoctial.f,
            orbit_target_equinoctial.g,
            orbit_target_equinoctial.h,
            orbit_target_equinoctial.k,
            orbit_target_equinoctial.L]

    # define histories to save outcomes to
    Q_hist = []
    orbit_hist = []
    time_hist = []
    mass_hist = []
    #initialize_values
    step_count = 0
    Q_value = Q(orbit, orbit_target, Q_Params, Sat_Params)
    time = time_initial
    mass = Sat_Params.mass
    # add initial values to history
    push!(Q_hist, Q_value)
    push!(orbit_hist, orbit)
    push!(time_hist, time)
    push!(mass_hist, mass)
    # integrate until hits max number of steps or Q convergence
    while (step_count < Q_Params.max_iter) && (Q_value > Q_Params.Q_convergence)
        thrust, F_t, F_r, F_n = evaluate_orbit_location(orbit, orbit_target, Q_Params, Sat_Params)
        if thrust
            GVE(orbit, time) = GaussVariationalEquationsEquinoctial(orbit, time; accel=[F_r,F_t,F_n])
            orbit_new, error = RK45(orbit, time, Q_Params.step_size, dOEdt = GVE)
            Q_new = Q(orbit_new, orbit_target, Q_Params, Sat_Params)
            if Q_new < Q_value
                Q_value = Q_new
                orbit = orbit_new
                mass = mass_loss(mass, Sat_Params, Q_Params.step_size)
            else
                orbit, error = RK45(orbit, time, Q_Params.step_size, dOEdt = GaussVariationalEquationsEquinoctial)
                Q_value = Q(orbit, orbit_target, Q_Params, Sat_Params)
            end
        else
            orbit, error = RK45(orbit, time, Q_Params.step_size, dOEdt = GaussVariationalEquationsEquinoctial)
            Q_value = Q(orbit, orbit_target, Q_Params, Sat_Params)
        end
        time = time + (Q_Params.step_size / 60. / 60. / 24.)
        push!(Q_hist, Q_value)
        push!(orbit_hist, orbit)
        push!(time_hist, time)
        push!(mass_hist, mass)
    end

    return(Q_hist, orbit_hist, time_hist, mass_hist)
end
