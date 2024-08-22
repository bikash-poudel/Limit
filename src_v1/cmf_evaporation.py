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
layer.wetness = 1

ET_pot = -1.005e-5 / 1000 / 1 * 1000 * 86400   #mm/day

ETpot = cmf.timeseries.from_scalar(ET_pot)  # constant ETpot of 20 mm/day (I know that this is too much)
et_pot_connection = cmf.timeseriesETpot(layer, cell.transpiration, ETpot)

# A solver, any is fine, really
solver = cmf.HeunIntegrator(p)
solver.t = datetime(2018, 5, 1)
et_act = cmf.timeseries(solver.t, cmf.day)
theta = cmf.timeseries(solver.t, cmf.day)
while solver.t < datetime(2018, 10, 1):
    et_act.add(cell.transpiration(solver.t))
    theta.add(layer.theta)
    solver(cmf.day)






