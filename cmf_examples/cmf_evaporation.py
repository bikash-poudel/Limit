import cmf
from datetime import datetime, timedelta

import matplotlib.pyplot as plt
# Plot results
from pylab import *

project = cmf.project()
# Add one cell at position (0,0,0), Area=1000m2 with a surface water storage
cell = project.NewCell(x=0, y=0, z=0, area=1000, with_surfacewater=False)

#r_curve = cmf.VanGenuchtenMualem(Ksat=0.106272, phi=0.35, alpha=0.193, n=2.0, m=0.099, theta_r=0.01)
# r_curve = cmf.VanGenuchtenMualem(Ksat=0.0106272, phi=0.35, alpha=0.005, n=2.22, m=0.25, theta_r=0.01)
# r_curve.w0 = 0.9995
#r_curve.l = 0.67
#r_curve.K = r_curve.Ksat * pow(r_curve.w0, 9.14)

r_curve = cmf.BrooksCoreyRetentionCurve(ksat=0.0106272, porosity=0.3500001, _b=9.14, theta_x=0.35, psi_x=-0.193)
r_curve.w0 = 0.9


L_boundaries = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 1.0, 1.25, 1.5,
                1.9, 2.3, 2.8, 3.5, 4.2, 5.0]

for i in range(20):
    depth = (i+1) * 0.05
    cell.add_layer(depth, r_curve)

cell.install_connection(cmf.Richards)

# weather
#summer = cmf.Weather(Tmax=30, Tmin=30, rH=20, wind=2, Rs=0)  # Rs==0 for 0 transpiration
#cell.set_weather(summer)

# Evaporation

# cell.install_connection(cmf.ShuttleworthWallace)
ETpot = cmf.timeseries.from_scalar(2.5)
cmf.timeseriesETpot(cell.layers[0], cell.evaporation, ETpot)
cell.vegetation.Height = 0
cell.vegetation.LAI = 0
cell.vegetation.CanopyClosure = 0


# Create the lower boundary condition
#gw = project.NewOutlet('groundwater',x=0,y=0,z=-2)
# Set the potential
#gw.potential = 0
# Connect the lowest layer to the groundwater using Richards percolation
#gw_flux=cmf.Richards(cell.layers[-1], gw)

# solver
solver = cmf.CVodeKrylov(project, 1e-6)
solver.t = cmf.Time(1, 1, 2011)

# Set all layers to a potential of -2 m
cell.saturated_depth = 0
# 100 mm water in the surface water storage for percolation
#cell.surfacewater.depth = 0.0

# Save potential and soil moisture for each layer,
# start with initial conditions
potential = [cell.layers.potential]
moisture = [cell.layers.theta]
ev_act, tr_act = [], []
flux = []
ql_L1 = []
# The run time loop, run for 7 days
for t in solver.run(solver.t,
                    solver.t + timedelta(days=250),
                    timedelta(days=1)):
    potential.append(cell.layers.potential)
    moisture.append(cell.layers.theta)
    ev_act.append(cell.evaporation(t))
    ql_L1.append(cell.layers[0].flux_to(cell.layers[1], t))
    tr_act.append(cell.transpiration(t))
    ql = [l.flux_to(next_l, t) for l, next_l in zip(cell.layers[:-1], cell.layers[1:])]

    flux.append(ql)
days = [50, 100, 150, 200, 250]
depth = - np.cumsum(cell.layers.thickness)
"""
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
"""

plot(ql_L1, label='ql_L1')
xlabel(r'$time [h]$')
ylabel(r'flux $[mm/day]$')
legend()
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

plot(moisture[50], depth, label='50 days')
plot(moisture[100], depth, label='100 days')
plot(moisture[150], depth, label='150 days')
plot(moisture[200], depth, label='200 days')
plot(moisture[250], depth, label='250 days')
xlabel(r'Soil moisture $\theta [m^3/m^3]$')
ylabel(r'$depth [cm]$')
plt.legend()
grid()
show()

plot(flux[49], depth[1:], label='50 days')
plot(flux[99], depth[1:], label='100 days')
plot(flux[149], depth[1:], label='150 days')
plot(flux[199], depth[1:], label='200 days')
plot(flux[249], depth[1:], label='250 days')
xlabel(r'moisture flux $[mm/day]$')
ylabel(r'$depth [cm]$')
plt.legend()
grid()
show()

th1 = [m[0] for m in moisture]
th2 = [m[1] for m in moisture]
th3 = [m[2] for m in moisture]
plot(th1, label='theta l1')
plot(th2, label='theta l2')
plot(th3, label='theta l3')
xlabel(r'$time [h]$')
ylabel(r'$\theta [m^3/m^3]$')
legend()
grid()
show()
