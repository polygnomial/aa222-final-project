using StaticArrays

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
    Ω # raan [rad]
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
    Ω = atan(k,h)
    ν = L - atan(g,f)
    return(KeplarianOrbit(a, e, i, ω, Ω, ν))
end

function Keplarian2Equinoctial(orbit::KeplarianOrbit)
    # unload orbit parameters
    a = orbit.a
    e = orbit.e
    i = orbit.i
    ω = orbit.ω
    Ω = orbit.Ω
    ν = orbit.ν
    # convert orbital elements
    p = a*(1-e^2)
    f = e*cos(ω+Ω)
    g = e*sin(ω+Ω)
    h = tan(i/2)*cos(Ω)
    k = tan(i/2)*sin(Ω)
    L = Ω + ω + ν
    return(EquinoctialOrbit(p,f,g,h,k,L))
end

function GaussVariationalEquationsEquinoctial(orbit, accel, time)
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

function RK45(orbit, GVE, step;dOEdt = GaussVariationalEquationsEquinoctial)
    k1 = step *
    k2 = step *
    k3 = step *
    k4 = step *
    k5 = step *
    k6 = step *
end
