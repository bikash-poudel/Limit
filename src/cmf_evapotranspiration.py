#loads libraries
import cmf
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
# from pylab import plot, twinx, ylabel, xlabel

# Create project with 1 cell and 1 water storage capacity
Braud_Ksat = 0.0106272  # in m/d
Braud_phi = 0.35  # porosity
Braud_alpha = 0.193  # Scale value of the water pressure(m)
Braud_n = 2.22
Braud_m = 0.099
Braud_theta_r = 0.01
Braud_eta = 9.14

VGM_retention_curve = cmf.VanGenuchtenMualem(Ksat=Braud_Ksat,
                                             phi=Braud_phi,
                                             alpha=Braud_alpha,
                                             n=Braud_n,
                                             m=Braud_m,
                                             theta_r=Braud_theta_r)

#Create a project
p = cmf.project()

# Create a cell at position (0,0,0) with 1000m2 size (making conversion from m3 to mm trivial)
cell = p.NewCell(0, 0, 0, 1000)
layer = cell.add_layer(0.1, VGM_retention_curve)

# Set summertime, when the living is easy... (1)
summer = cmf.Weather(Tmin=19, Tmax=28, rH=50, wind=4.0,
                     sunshine=0.9, daylength=14, Rs=26)
cell.set_weather(summer)

# Initial condition (4)
layer.volume = 350 * layer.thickness # theta
# ET-Method (3)
et_pot_turc = cmf.TurcET(layer, cell.transpiration)



# Stress conditions (5, 6)
stress = cmf.VolumeStress(300, 100)
cell.set_uptakestress(stress)
# A solver, any is fine, really
solver = cmf.HeunIntegrator(p)
solver.t = datetime(2018, 5, 1)
et_act = cmf.timeseries(solver.t, cmf.day)
volume = cmf.timeseries(solver.t, cmf.day)
theta = cmf.timeseries(solver.t, cmf.day)
while solver.t < datetime(2018, 10, 1):
    et_act.add(cell.transpiration(solver.t))
    volume.add(layer.volume)
    theta.add(layer.theta)
    solver(cmf.day)

# And a plot
plt.plot(et_act, c='g')
plt.ylabel(r'$ET_{act} \left[\frac{mm}{day}\right]$')
plt.twinx()
plt.plot(theta, c='b')
plt.ylabel('Theta in mm')
plt.show()
