# RKF45 constants (to be stored in seperate file)
const A2 = 2. / 9.
const A3 = 1. / 3.
const A4 = 3. / 4.
const A5 = 1.
const A6 = 5. / 6.
const B21 = 2. / 9.
const B31 = 1. / 12.
const B41 = 69. / 128.
const B51 = -17. / 12.
const B61 = 65. / 432.
const B32 = 1. / 4.
const B42 = -243. / 128.
const B52 = 27. / 4.
const B62 = -5. / 16.
const B43 = 135. / 64.
const B53 = -27. / 5.
const B63 = 13. / 16.
const B54 = 16. / 15.
const B64 = 4. / 27.
const B65  = 5. / 144.
const CH1 = 47. / 450.
const CH2 = 0.
const CH3 = 12. / 25.
const CH4 = 32. / 225.
const CH5 = 1. / 30.
const CH6 = 6. / 25.
const CT1 = -1. / 150.
const CT2 = 0.
const CT3 = 3. / 100.
const CT4 = -16. / 75.
const CT5 = -1. / 20.
const CT6 = 6. / 25.

# https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
function RK45(orbit, time, step;dOEdt = GaussVariationalEquationsEquinoctial)
    # compute k1 through k6
    k1 = step * dOEdt(orbit, time)
    modifier = B21 * k1
    k2 = step * dOEdt(orbit + modifier, time + (A2 * step))
    modifier = B31 * k1 + B32 * k2
    k3 = step * dOEdt(orbit + modifier, time + (A3 * step))
    modifier = B41 * k1 + B42 * k2 + B43 * k3
    k4 = step * dOEdt(orbit + modifier, time + (A4 * step))
    modifier = B51 * k1 + B52 * k2 + B53 * k3 + B54 * k4
    k5 = step * dOEdt(orbit + modifier, time + (A5 * step))
    modifier = B61 * k1 + B62 * k2 + B63 * k3 + B64 * k4 + B65 * k5
    k6 = step * dOEdt(orbit + modifier, time + (A6 * step))
    # compute weighted average
    orbit_new = orbit + CH1 * k1 + CH2 * k2 + CH3 * k3 + CH4 * k4 + CH5 * k5 + CH6 * k6
    # compute error
    error = abs.(CT1 * k1 + CT2 * k2 + CT3 * k3 + CT4 * k4 + CT5 * k5 + CT6 * k6)
    # NOTE: we can add in adaptive time steps here if we want...
    # return new orbit values and numerical integration error
    return(orbit_new, error)
end
