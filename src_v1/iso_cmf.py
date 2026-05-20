# loads libraries
import cmf
from datetime import datetime, timedelta
from math import exp
import numpy as np
from numpy.random import rand

from iso_time import iso_time
from iso_layers import iso_layer
from Boundary_conditions import BoundaryCondition, iso_atmosphere
import solve_iso_transport as iso
from Visualize import Visualize

import matplotlib.pyplot as plt

# Setup CMF
# A function to create a retention curve (used for the whole profile)
# for a specific depth, using an exponential decline function for Ksat with depth

# Likelihood Tests (Mathieu and Bariac (1996).)
ignorealphai = True
ignorealphaik = True
ignoredl = False
ignoredv = True

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

VGM_retention_curve.w0 = 0.9995  # See: https://philippkraft.github.io/cmf/cmf_tut_retentioncurve.html Oversaturation
# Create a project
p = cmf.project()

# Create a cell at position (0,0,0) with 1000m2 size (making conversion from m3 to mm trivial)
c = p.NewCell(0, 0, 0, 1000)
#c.vegetation.RootDepth = 0.0

# Customize cell
# Top layer thickness of e10-5 m as per SISPAT_iso
lower_boundaries_of_layer = [0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                             0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]


for d in lower_boundaries_of_layer:
    l = c.add_layer(d, VGM_retention_curve)
    #c.saturated_depth = d   # initial wetness gives wrong values of theta and wetness. (theta > porosity,  wetness larget then 1)
    l.wetness = 1   ##### Check tehts exceeds porosity if  l.wetness = 1.0
    #c.saturated_depth = 0  # check for wetness = 1, saturated depth = 0,

# Create a surfacewater storage
#c.surfacewater_as_storage()

# use Richards connection
c.install_connection(cmf.Richards)



#Evaportation
ET_CMF = True
E_pot = -1.005e-5  # kg/(m2 * s) or[ mm s-1]

if ET_CMF == False:
    # Create a Neumann Boundary condition connected to top layer to account for flows due to evaporation
    E_act = cmf.NeumannBoundary_create(c.layers[0])
    E_act.flux = cmf.timeseries_from_scalar(0.0)
elif ET_CMF == True: #use CMF to calculate the ETact based on a fixed ETpot
    potET = abs(E_pot) * 86400  # mm/day
    ETpot = cmf.timeseries.from_scalar(potET)

    for layer in c.layers:
        cmf.timeseriesETpot(layer, c.evaporation, ETpot)  # connection to top only or to all layers ???
        # et_pot_connection = cmf.timeseriesETpot(layer, c.evaporation, ETpot)
else:
    raise NotImplementedError

# Simulation period
start = cmf.Time(1, 1, 2020)
end = start + timedelta(days=50)  # run cmf for 300 days
dt = cmf.h  # time step

# Make solver
# Different solvers in cmf2
solver = cmf.CVodeBanded(p, 1e-09)
#solver.t = start
# Complete Jacobian, alternatives: CVodeBanded, for komplex sparse Jacobians: CVodeKLU,
# CVodeKrylov -> like Banded but with "Krylov preconditioner", CVodeAdams an very advance explicit solver

solver.t = start

# Setup
initial_c_2H_soil = iso.delta_to_concentration(delta_i=-65, solute_i="2H")
initial_c_18O_soil = iso.delta_to_concentration(delta_i=-8, solute_i="18O")
initial_c_2H_atmosphere = iso.delta_to_concentration(delta_i=-65, solute_i="2H")
initial_c_18O_atmosphere = iso.delta_to_concentration(delta_i=-8, solute_i="18O")

T_atmosphere = 303.0
rH_atmosphere = 0.20

my_iso_atmosphere = iso_atmosphere(
    initial_c_atmosphere={"2H": initial_c_2H_atmosphere, "18O": initial_c_18O_atmosphere},
    # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
    initial_T_atmosphere=T_atmosphere,  # temperature in the atmosphere in Kelvin
    initial_Rh_atmosphere=rH_atmosphere)  # relative humidity of the atmosphere (-)

# my_iso_cell = iso_cell(atmosphere=my_iso_atmosphere, area=1.0)

i_upper_boundary = 0.0
my_layers = []
for i_lower_boundary, i_cmf_layer in zip(lower_boundaries_of_layer, c.layers):
    my_iso_layer = iso_layer(upper_boundary=i_upper_boundary,
                             lower_boundary=i_lower_boundary,
                             initial_c_solutes={"2H": initial_c_2H_soil, "18O": initial_c_18O_soil},
                             # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                             initial_theta=i_cmf_layer.theta,  # volumetric water content m3/m3
                             theta_0=0.01,  # volumetric water content m3/m3 at high suctions
                             theta_sat=i_cmf_layer.porosity,  # volumetric water content m3/m3 at saturation
                             initial_T=303,  # soil temperature in Kelvin
                             porosity=i_cmf_layer.porosity,  # porosity of the soil m3/m3
                             tortuosity=0.67)  # tortuosity of the soil m/m
    i_upper_boundary = i_lower_boundary
    my_layers.append(my_iso_layer)

my_time = iso_time(final_time=(end-start)/cmf.h,
                   delta_time=1,
                   time_units='hours')


# The run time loop. Iterates over the outer timestep of the model
Ciso = {'2H': [], '18O': []}
#et_act = cmf.timeseries(solver.t, cmf.day)

theta_t = [np.array(c.layers.theta)]
potential = [c.layers.potential]
moisture = [c.layers.theta]
hs = []
ql_up_t = []
ql_down_t = []
evap = []
while solver.t < end:

    # run cmf for one time step
    # solver.reset() # - Resets the solver history. Only needed if you really change something in the system
    solver(solver.t + dt)

    theta_t1 = c.layers.theta
    theta_t.append(theta_t1)

    ql_layers_up_t1 = []
    ql_layers_down_t1 = []
    qv_up = []
    qv_down = []

    f = cmf.sec/cmf.day / c.area  # convert to m per second
    qqq = 0.0
    # Todo: Implement like this:
    for i_iso_layer, i_cmf_layer in zip(my_layers, c.layers):
        if i_cmf_layer.upper is None and i_cmf_layer.lower is not None:
            ql_layers_up_t1.append(0.0 * f + qqq)
            ql_layers_down_t1.append(i_cmf_layer.flux_to(i_cmf_layer.lower, solver.t) * f + qqq)
        elif i_cmf_layer.upper is not None and i_cmf_layer.lower is not None:
            ql_layers_up_t1.append(i_cmf_layer.flux_to(i_cmf_layer.upper, solver.t) * f + qqq)
            ql_layers_down_t1.append(i_cmf_layer.flux_to(i_cmf_layer.lower, solver.t) * f + qqq)
        elif i_cmf_layer.upper is not None and i_cmf_layer.lower is None:
            ql_layers_up_t1.append(i_cmf_layer.flux_to(i_cmf_layer.upper, solver.t) * f + qqq)
            ql_layers_down_t1.append(0.0 * f)

    qv_up = [0] * len(my_layers)
    qv_down = [0] * len(my_layers)

    ql_up_t.append(ql_layers_up_t1)
    ql_down_t.append(ql_layers_down_t1)

    def soil_relative_humidity(h, T, R=461.5, g=9.81):

        """
        Calculates the relative humidity at given soil depth for given pressure head and temperature  certain
        soil depth according to Braud et al. 2005 Eq. 46

        Description
        == == == == == =

         @param h: pressure head [m]

         @param R: perfect gas mass constant for water vapour [J kg-1]

         @param T: Temperaturtre in kelvin [k]

         @param g: acceleration due to gravity [ms-2]
        """

        hu = exp(g * h / (R * T))

        return hu

    # Evaporation

    # calculate actual Evaporation for the next time step
    ha_surface = rH_atmosphere * iso_atmosphere.pv_sat(T_atmosphere) / iso_atmosphere.pv_sat(my_layers[0].T)
    hs_surface = soil_relative_humidity(c.layers[0].matrix_potential, 303)

    if ET_CMF == False:
        actual_Evaporation = E_pot / 1000 * c.area * 86400 * (hs_surface - ha_surface) / (1 - ha_surface)   # m3 per day
    elif ET_CMF == True:
        actual_Evaporation = c.evaporation(solver.t)

    evap.append(actual_Evaporation)
    hs.append(hs_surface)

    potential.append(c.layers.potential)
    moisture.append(c.layers.theta)

    fluxes = [ql_layers_up_t1, ql_layers_down_t1, qv_up, qv_down, [[0.35]*len(my_layers), [0.35]*len(my_layers)]]

    def BC(evap):
        # Boundary conditions
        U_boundary_conc = iso.delta_to_concentration(delta_i=-70, solute_i='2H')
        L_boundary_conc = iso.delta_to_concentration(delta_i=-65, solute_i='2H')

        BC = BoundaryCondition()
        BC.upper_boundary('dirichlet', U_boundary_conc)  # neuman[flux] / dirichlet[constant conc] / atmosphere[potential evaporation]
        BC.lower_boundary('neuman', 0)  # TODO: need to check lower boundary

        return BC


    B_C = BC(E_pot)

    delta_t = 3600  # descrete time dt [s]
    time_steps = my_time.total_seconds() / delta_t

    for t_steps in range(int(time_steps)):
        # Run iso_top simulation
        Ci_t = iso.run_1D_model(my_iso_atmosphere, my_layers, fluxes, B_C, my_time, hs_surface, solutes=['2H'],
                                ignore_alpha_i=ignorealphai,
                                ignore_alpha_i_k=ignorealphaik,
                                ignore_dl_i=ignoredl,
                                ignore_dv_i=ignoredv)

        update_layers = iso.update_c_i(Ci_t['2H'], '2H', my_layers)
        my_layers = update_layers
        # my_layers = iso.update_c_i(Ci_t['2H'], '18O', my_layers)

    Ciso['2H'].append(Ci_t['2H'])
    # Ciso['18O'].append(Ci_t['18O'])


C2H = np.array(Ciso['2H']).T
#C18O = np.array(Ciso['18O']).T

theta = np.array(theta_t).T
ql_up = np.array(ql_up_t).T
ql_down = np.array(ql_down_t).T

# Convert Result into delta notation
my_vectorized_function = np.vectorize(iso.concentration_to_delta, excluded=['solute_i'])
C_delta = my_vectorized_function(C2H, solute_i='2H')
#C18O_delta = my_vectorized_function(C18O, solute_i='18O')


#Visualize
plot = Visualize(my_layers, my_time)
plot.profile(C_delta, solute_i='2H', print_time_steps=24*5)
plot.breakthrough(C_delta, solute_i='2H', print_steps=2)

plt.plot(evap, label='evaporation')
plt.legend()
plt.show()

plt.plot(ql_up[:, -1], -np.arange(len(my_layers)), label='flux_up at final time')
#plt.plot(ql_down[:, -1], -np.arange(len(my_layers)), label='flux_down at final time')
plt.legend()
plt.show()

plt.plot(theta[:, -1], -np.arange(len(my_layers)), label='soil moisture at the end')
plt.legend()
plt.show()

plt.plot(theta[0], label='theta_top_layer')
plt.legend()
plt.show()

plt.plot(hs, label='soil surface relative humidity')
plt.legend()
plt.show()

plt.plot(potential)
plt.title('potential')
plt.show()

plt.plot(moisture)
plt.title('moisture')
plt.show()

