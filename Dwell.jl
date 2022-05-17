include("Qlaw.jl")
include("Integrator.jl")
using Plotly

# Create starting orbit and time, end time, step size
orbit_inital_keplarian = KeplarianOrbit(
                450. + R_e,
                0.0001,
                deg2rad(0.001),
                deg2rad(0.001),
                deg2rad(0.001),
                deg2rad(0)
                )
orbit_initial_equinoctial = Keplarian2Equinoctial(orbit_inital_keplarian)
orbit_initial = [orbit_initial_equinoctial.p,
                orbit_initial_equinoctial.f,
                orbit_initial_equinoctial.g,
                orbit_initial_equinoctial.h,
                orbit_initial_equinoctial.k,
                orbit_initial_equinoctial.L]

time_initial = 59715.0 # MJD
time_end = time_initial + 0.1 # stop integrator after this time
step_size = (1. / 24.) / 60. / 60. # 1 min

# Initialize history to store propogated orbits and error bounds
orbit_history = [orbit_initial]
error_history = [[0.,0.,0.,0.,0.,0.]]
time_history = [time_initial]

# Integrate orbit with no external forces
time = time_initial
print("Starting integration...\n")
while time < time_end
    global time = time_history[end]
    orbit = orbit_history[end]
    orbit_new, error = RK45(orbit, time, step_size, dOEdt = GaussVariationalEquationsEquinoctial)
    push!(time_history, time + step_size)
    push!(orbit_history, orbit_new)
    push!(error_history, error)
end
print("Integration finished!\n")
# store ECI terms
ECI_x = []
ECI_y = []
ECI_z = []
for orbit in orbit_history
    equinoctial_orbit = EquinoctialOrbit(orbit[1], orbit[2], orbit[3], orbit[4], orbit[5], orbit[6])
    ECI_orbit = Equinoctial2ECI(equinoctial_orbit)
    push!(ECI_x, ECI_orbit[1])
    push!(ECI_y, ECI_orbit[2])
    push!(ECI_z, ECI_orbit[3])
end

# plot orbit
n = 10
u = collect(range(0,2*π, length = n))
v = collect(range(0,π, length = n))
r_e = 6378.

u_range = range(0, stop=2π, length=100)
u = u_range' .* ones(100)
v_range = range(0, stop=π, length=100)
v = ones(100)' .* v_range

x = @. r_e * cos(u) * sin(v)
y = @. r_e * sin(u) * sin(v)
z = @. r_e * cos(v)

layout = Layout(title="Demonatration of Numerical Propogator")
trace1 = surface(x=x, y=y, z=z, showscale=false)
trace2 = scatter3d(x=(ECI_x ./ 1000), y=(ECI_y ./ 1000), z=(ECI_z ./ 1000),
    mode="markers",
    marker=attr(
        size=1,
    ))
p = plot([trace1, trace2], layout)
savefig(p, "orbit_propogator.png")
# Propogate orbit with known solutions to determine error
#TBD!!!
# Plot differences and compare to error from integrator
#TBD!!!
