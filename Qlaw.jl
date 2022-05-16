using StaticArrays
using LinearAlgebra

# constants
µ = 3.986004419 * 10^14

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

function Equinoctial2Keplarian(orbit::KeplarianOrbit)
    # unload orbit parameters
    p = orbit.p
    f = orbit.f
    g = orbit.g
    h = orbit.h
    k = orbit.k
    L = orbit.L
    # convert orbital elements
    # from https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
    a = p/(1-f^2-g^2)
    e = sqrt(f^2 + g^2)
    i = atan(2*sqrt(h^2+k^2), 1-h^2-k^2)
    ω = atan(g*h-f*k, f*h+g*k)
    Ω = atan(k,h)
    ν = L - atan(g,f)
    return(KeplarianOrbit(a, e, i, ω, Ω, ν))
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

function GaussVariationalEquationsEquinoctial(orbit, accel)
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
    A = [[0, sinL, -cosL, 0, 0, 0]
         [2*p/w, (1/w)*((w+1)*cosL+f), (w+1)*sinL+g, 0, 0, 0]
         [0, -(g/w)*(h*sinL-k*cosL), (f/w)*(h*sinL-k*cosL), (s_2*cosL)/(2*w), (s_2*sinL)/(2*w), h*sinL-k*cosL]]
    A = []
    A = root_p_µ * A
    b = [0, 0, 0, 0, 0, sqrt(µ*p)*(w/p)^2]

    # calculate doe/dt
    dOEdt = A
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
    print(Q_params)
    # unpack orbit
    p = orbit[1] # semi-latus rectum [meters]
    f = orbit[2]
    g = orbit[3]
    h = orbit[4]
    k = orbit[5]
    L = orbit[6] # true longitude [rad]
    s2 = 1+h^2+k^2 # derived
    
    p_t = target[1] # semi-latus rectum [meters]
    f_t = target[2]
    g_t = target[3]
    h_t = target[4]
    k_t = target[5]
    L_t = target[6] # true longitude [rad]
    
    # unpack params
    r_p_min = Q_params.r_p_min
    W = zeros(5)
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

    a = p/(1-f^2-g^2) # semi-major axis
    a_t = p_t/(1-f_t^2-g_t^2)
    e = sqrt(f^2 + g^2) # eccentricity
    r_p = p/(1+e) # periapsis

    #         
    sqrt_fg = sqrt(f^2 + g^2)
    sqrt_p_over_µ = sqrt(p/µ)
    oe_dot_xx_temp = [
        2*a*sqrt((a*(1 + sqrt_fg))/(µ*(1 - sqrt_fg))),
        2*sqrt_p_over_µ,
        0.0,
        1/2*sqrt_p_over_µ*s2/(sqrt(1-g^2)+f),
        1/2*sqrt_p_over_µ*s2/(sqrt(1-f^2)+g)
    ]
    oe_dot_xx_temp[3] = oe_dot_xx_temp[2]
    oe_dot_xx = oe_dot_xx_temp * F # need to 
    
    doe_dot_xx_dF = [0.0, 0.0, 0.0]
    # avoid divide by zero if F is aligned with a coordinate axis
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
    
    # Scaling Factors
    S_oe = [(1 + ((1 - a_t)/(m*a_t))^n)^(1/r), # eq. 8
        1,
        1,
        1,
        1
    ]

    # Penalty
    P = exp(k_penalty*(1 - r_p/r_p_min))
        
    Q_p = zeros(5)
    Q = zeros(5)
    Q_m = zeros(5)
    
    Q_p[1] = S_oe[1]*W[1]*((a + h/2 - target[1])/oe_dot_xx[1])^2
    Q_p[2] = S_oe[2]*W[2]*((orbit[2] + h/2 - target[2])/oe_dot_xx[2])^2
    Q_p[3] = S_oe[3]*W[3]*((orbit[3] + h/2 - target[3])/oe_dot_xx[2])^2
    Q_p[4] = S_oe[4]*W[4]*((orbit[4] + h/2 - target[4])/oe_dot_xx[3])^2
    Q_p[5] = S_oe[5]*W[5]*((orbit[5] + h/2 - target[5])/oe_dot_xx[4])^2
        
    Q[1] = S_oe[1]*W[1]*((a - target[1])/oe_dot_xx[1])^2
    Q[2] = S_oe[2]*W[2]*((orbit[2] - target[2])/oe_dot_xx[2])^2
    Q[3] = S_oe[3]*W[3]*((orbit[3] - target[3])/oe_dot_xx[2])^2
    Q[4] = S_oe[4]*W[4]*((orbit[4] - target[4])/oe_dot_xx[3])^2
    Q[5] = S_oe[5]*W[5]*((orbit[5] - target[5])/oe_dot_xx[4])^2
    
    Q_m[1] = S_oe[1]*W[1]*((a - h/2 - target[1])/oe_dot_xx[1])^2
    Q_m[2] = S_oe[2]*W[2]*((orbit[2] - h/2 - target[2])/oe_dot_xx[2])^2
    Q_m[3] = S_oe[3]*W[3]*((orbit[3] - h/2 - target[3])/oe_dot_xx[2])^2
    Q_m[4] = S_oe[4]*W[4]*((orbit[4] - h/2 - target[4])/oe_dot_xx[3])^2
    Q_m[5] = S_oe[5]*W[5]*((orbit[5] - h/2 - target[5])/oe_dot_xx[4])^2
        
    dQ = sum(4 .* Q .+ (1 + W[1] .* P)) * sum((Q_p - Q_m)/(h)) # TODO double check this sum is correct
    
    return dQ .* doe_dot_xx_dF

end