using StaticArrays
using LinearAlgebra

# constants
const µ = 3.986004419 * 10^14 # standard gravitational parameter earth [m^3/s^2]
const R_e = 6371000. # average earth radius [m]

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
function Equinoctial2Keplarian(orbit::KeplarianOrbit)
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
    r_p_min
    W_a # weights
    W_f
    W_g
    W_h
    W_k
    m # penalty params
    n
    r
    k # I forget
    h # central difference stepsize
end

function calc_D(orbit, target, F_vec, µ, Q_params)
    # unpack orbit
    a = orbit[1]
    f = orbit[2]
    g = orbit[3]
    h = orbit[4]
    k = orbit[5]

    # derived elements
    s2 = 1+h^2+k^2
    e = sqrt(f^2 + g^2) # eccentricity
    p = a*(1-e^2) # needed for derivatives (eq. 14-18)
    r_p = p/(1+e) # periapsis

    # hyper params
    r_p_min = Q_params.r_p_min
    W = @MVector zeros(5)
    W[1] = Q_params.W_a
    W[2] = Q_params.W_f
    W[3] = Q_params.W_g
    W[4] = Q_params.W_h
    W[5] = Q_params.W_k
    m = Q_params.m # penalty params
    n = Q_params.n
    r = Q_params.r
    k_penalty = Q_params.k # I forget
    h = Q_params.h # central difference stepsize

    # F
    F = norm(F_vec)

    # Eq. 14-18
    sqrt_fg = sqrt(f^2 + g^2) # helper
    sqrt_p_over_µ = sqrt(p/µ) # helper
    oe_dot_xx_temp = [
        2*a*sqrt((a*(1 + sqrt_fg))/(µ*(1 - sqrt_fg))),
        2*sqrt_p_over_µ,
        0.0,
        1/2*sqrt_p_over_µ*s2/(sqrt(1-g^2)+f),
        1/2*sqrt_p_over_µ*s2/(sqrt(1-f^2)+g)
    ]
    oe_dot_xx_temp[3] = oe_dot_xx_temp[2]
    oe_dot_xx = oe_dot_xx_temp * F # need to


    # Eq. 21-23 - doe_dot/dF (3 partials)
    # avoid divide by zero if F is aligned with a coordinate axis
    doe_dot_xx_dF = @MVector zeros(3)
    if F_vec[2] == 0 && F_vec[3] == 0
        doe_dot_xx_dF[1] = sum(oe_dot_xx_temp) * F_vec[1]
    else
        doe_dot_xx_dF[1] = sum(oe_dot_xx_temp) * F_vec[1] / F
    end
    if F_vec[1] == 0 && F_vec[3] == 0
        doe_dot_xx_dF[2] = sum(oe_dot_xx_temp) * F_vec[2]
    else
        doe_dot_xx_dF[2] = sum(oe_dot_xx_temp) * F_vec[2] / F
    end
    if F_vec[1] == 0 && F_vec[2] == 0
        doe_dot_xx_dF[3] = sum(oe_dot_xx_temp) * F_vec[3]
    else
        doe_dot_xx_dF[3] = sum(oe_dot_xx_temp) * F_vec[3] / F
    end

    # Eq. 8 - Scaling Factors
    a_t = target[1]
    S_oe = @MVector [(1 + ((1 - a_t)/(m*a_t))^n)^(1/r), # eq. 8
        1,
        1,
        1,
        1
    ]

    # Eq. 9 - Penalty
    P = exp(k_penalty*(1 - r_p/r_p_min))

    Q_p = @MVector zeros(5)
    Q_m = @MVector zeros(5)

    # Eq. 7 - dQ/doe
    Q_p = S_oe .* W .*((orbit .+ h/2 .- target)./oe_dot_xx).^2
    Q_m = S_oe .* W .*((orbit .- h/2 .- target)./oe_dot_xx).^2
    dQ = sum(1 + W[1] .* P) * sum((Q_p .- Q_m)/h) # TODO double check this sum is correct

    return dQ .* doe_dot_xx_dF

end
