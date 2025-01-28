import cmf
from datetime import datetime,timedelta

# Plot results
from pylab import *

project = cmf.project()
# Add one cell at position (0,0,0), Area=1000m2 with a surface water storage
cell = project.NewCell(x=0,y=0,z=0,area=1000, with_surfacewater=True)

"""
r_curve = cmf.VanGenuchtenMualem(Ksat=0.0106272, phi=0.35, alpha=0.01, n=2.0)
r_curve.w0 = 0.9995
r_curve.l = 0.67
r_curve.K = r_curve.Ksat * pow(r_curve.w0, 9.14)
"""

r_curve = cmf.BrooksCoreyRetentionCurve(ksat=0.0106272, porosity=0.35, _b=1/0.22, theta_x=0.3499)
r_curve.w0 = 0.9995

for i in range(10):
    depth = (i+1) * 0.1
    cell.add_layer(depth, r_curve)

cell.install_connection(cmf.Richards)

# Create the boundary condition
gw = project.NewOutlet('groundwater',x=0,y=0,z=-2)
# Set the potential
gw.potential = 0
# Connect the lowest layer to the groundwater using Richards percolation
gw_flux=cmf.Richards(cell.layers[-1],gw)

solver = cmf.CVodeKrylov(project,1e-6)
solver.t = cmf.Time(1,1,2011)

# Set all layers to a potential of -2 m
cell.saturated_depth = 0
# 100 mm water in the surface water storage for percolation
cell.surfacewater.depth = 0.1

# Save potential and soil moisture for each layer,
# start with initial conditions
potential = [cell.layers.potential]
moisture = [cell.layers.theta]
ev_act, tr_act = [], []
flux = []
# The run time loop, run for 7 days
for t in solver.run(solver.t,
                    solver.t + timedelta(days=7),
                    timedelta(hours=1)):
    potential.append(cell.layers.potential)
    moisture.append(cell.layers.theta)
    ev_act.append(cell.evaporation(t))
    tr_act.append(cell.transpiration(t))
    ql = [l.flux_to(next_l, t) for l, next_l in zip(cell.layers[:-1], cell.layers[1:])]
    flux.append(ql)

subplot(211)
plot(moisture)
ylabel(r'Soil moisture $\theta [m^3/m^3]$')
xlabel(r'$time [h]$')
grid()
subplot(212)
plot(potential)
ylabel(r'Water head $\Psi_{tot} [m]$')
xlabel(r'$time [h]$')
grid()
show()

plot(ev_act, label='evap')
plot(tr_act, label='trans')
plot(np.array(ev_act) + np.array(tr_act), label='qev + qtr')
xlabel(r'$time [h]$')
ylabel(r'flux $[mm/day]$')
legend()
grid()
show()

plot(moisture[-1], -np.arange(0, len(moisture[0])))
xlabel(r'Soil moisture $\theta [m^3/m^3]$')
ylabel(r'$depth [cm]$')
grid()
show()

plot(flux[-1], -np.arange(1, len(moisture[0])))
xlabel(r'moisture flux $[mm/day]$')
ylabel(r'$depth [cm]$')
grid()
show()

