
# -*- coding: utf-8 -*-
"""
one-dimensional numerical solution for liquid phase isotope transport

units: meter - second - kilogram
"""

import sys
import math
import os
import datetime
import numpy as np
import scipy
import scipy.special as ssp
import matplotlib.pyplot as plt
import pathlib
from scipy.optimize import fsolve

__author__ = "Bikash Poudel (bikash.poudel@outlook.com)"
__version__ = "$Revision: 0.01 $"
__date__ = datetime.date(2019, 9, 5)
__copyright__ = "Copyright (c) 2019 Bikash Poudel"
__license__ = "BSD3"


def space_time(th, a, v, D, c0, col_dx, col_nx, time_dt, n_timesteps, L_BC, R_BC):
    # initialisation

    space_time_array = np.zeros(shape=(n_timesteps, col_nx))
    # coefficient matrix
    coeff = np.zeros(shape=(col_nx, col_nx))
    f = np.zeros(shape=(col_nx))  # weighting function
    # coeff[0,:] = c0 # at t=0 some initial concentration
    # space_time_array[0,:] = 1 # at t=0 some initial concentration

    # weighting functions
    a1 = - th * ((a * v / col_dx) + (D / (col_dx * col_dx)))
    a2 = (1 / time_dt) + (2 * th * D / (col_dx * col_dx)) - (th * v * (1 - (2 * a)) / col_dx)
    a3 = th * (((1 - a) * v / col_dx) - (D / (col_dx * col_dx)))
    # print("a1=",a1,"a2=",a2,"a3=",a3)

    f1 = (1 - th) * ((D / (col_dx * col_dx)) + (a * v / col_dx))
    f2 = (1 / time_dt) - (2 * (1 - th) * D / (col_dx * col_dx)) + ((1 - th) * (1 - (2 * a)) * v / col_dx)
    f3 = (1 - th) * ((D / (col_dx * col_dx)) - ((1 - a) * v / col_dx))
    # print("f1 =",f1,"f2 = ",f2,"f3 = ",f3)

    if L_BC == 0:  # Left Boundry Condition
        coeff[0, 0] = 1
        f[0] = c0  # constant concentration
        # coeff[0,1] = a3
        """  
        else:
        coeff[0,0] = a2 - (2 * a1 * v * col_dx / D )
        coeff[0,1] = a1 + a2
        f[0] = 
        """

    if R_BC == 0:  # Right Boundary Condition
        coeff[col_nx - 1, col_nx - 2] = a1
        coeff[col_nx - 1, col_nx - 1] = a2 + a1  # Cn+1 = Cn

    for i in range(1, col_nx - 1):
        coeff[i, i - 1] = a1
        coeff[i, i] = a2
        coeff[i, i + 1] = a3
    # print("Coeff Matrix \n",coeff)
    for t in range(0, n_timesteps - 1):

        for s in range(1, col_nx - 1):  # new f matrix
            f[s] = space_time_array[t, s - 1] * f1 + space_time_array[t, s] * f2 + space_time_array[t, s + 1] * f3
        # if t == 1:
        # print("f \n ",f)

        C_t = np.linalg.solve(coeff, f)

        for i in range(0, col_nx - 1):
            space_time_array[t + 1, i] = C_t[i]  # Conc at t+1 th time in global space time matrix

        # print("f \n ",f)
    # print("space time array\n", space_time_array)

    return space_time_array


def main():
    soil_depth = 1.  # m
    n_layers = 10  # No of descrete layers
    layer_dx = soil_depth / n_layers  # m

    time_final = 60. * 60. * 24. * 20.  # days
    time_dt = 60. * 60. * 24.
    n_timesteps = int(time_final / time_dt)

    dh = soil_depth  # hydraulic head difference
    i = dh / soil_depth  # head gradient

    # Soil Hydraulic Parameters
    k_sat = 0.6  # conductivity m/day
    porosity = 0.6
    bulk_density = 1800  # kg/cum
    theta = 0.35  # for saturated flow
    tortuosity = 0.67

    hydrodynamic_dispersivity = 0

    q_l = k_sat * i  # flux
    #q = 0
    #v = (q / porosity) / (60. * 60. * 24.)  # physical velocity m per day
    #print("velocity = ", v * (60. * 60. * 24.))
    print("flux = ", q_l) #m/day

    def dl_i(tortuosity, theta,T=283.15, Isotopologue = '2H', formulation = "Cuntz", ignore_dl_i = False):

        try:

            if formulation == "Melayah":
                if "18O" == Isotopologue:
                    a_i = 0.9669
                elif "2H" == Isotopologue:
                    a_i = 0.9833
                else:
                    raise NotImplementedError

                dl_i = a_i * 1.10 ** -9 * math.exp(-(535400 / T ** 2) + (1393.3 / T) + 2.1876)

            # SLI: cable_sli_solve.f90::L2381-2390
            elif formulation == "Cuntz":
                if "18O" == Isotopologue:
                    a_i = 1 / 1.026
                elif "2H" == Isotopologue:
                    a_i = 1 / 1.013
                else:
                    raise NotImplementedError

                dl_i = a_i * 1e-7 * math.exp(-577 / (T - 145))

            else:
                raise NotImplementedError

            return dl_i * 86400 * theta * tortuosity * theta

        except ValueError as err:
            print(err)
            raise NotImplementedError

    dl_i_eff = dl_i(tortuosity, theta) + hydrodynamic_dispersivity * abs(q_l)
    print('effective liquid diffusion :', dl_i_eff)


    #alpha = 1
    #D = alpha * v + diffusion  # sq.m. per day Dispersion

    #k = 0  # decay coefficient

    # Boundary Conditions
    Upper_BC = 0  # 0 = constant concentration, 1 = Concentration flux
    Lower_BC = 0  # 0 = zero concentration gradient, 0.5 = constant concentration, 1 = Concentration Flux

    # initial condition
    c0 = 1  # initial concentration
    c_end = 0  # concentration at outlet

    # print information
    t_ime = 0.5  # of total time
    time = int(round(t_ime * (n_timesteps - 1)))

    p = 0.5  # of length
    pnt = int(round(p * (layer_dx - 2)))  # print at

    # Choose Time Descretization Scheme
    th = 0.5  # theta 0-explicit, 0.5-central, 1-implicit  in time

    # Numerical solutions using the weighting functions
    
    """
    space_time_array_forward = space_time(th, 0, q_l, dl_i_eff, c0, layer_dx, n_layers, time_dt, n_timesteps, Upper_BC, Lower_BC)
    space_time_array_central = space_time(th, 0.5, q_l, dl_i_eff, c0, layer_dx, n_layers, time_dt, n_timesteps, Upper_BC, Lower_BC)
    space_time_array_backward = space_time(th, 1, q_l, dl_i_eff, c0, layer_dx, n_layers, time_dt, n_timesteps, Upper_BC, Lower_BC)

    plt.subplot(211)
    plt.title("c/c0 vs. distance")
    plt.plot(space_time_array_forward[time, :],label='forward')
    plt.plot(space_time_array_central[time, :], label='central')
    plt.plot(space_time_array_backward[time, :], label='backward')
    plt.legend(loc='best')
    plt.xlabel("distance")
    plt.ylabel("concentration")
    plt.grid()
    # plt.show()

    plt.subplot(212)
    plt.title("c/c0 vs. time ")
    plt.plot(space_time_array_forward[:, pnt],label='forward')
    plt.plot(space_time_array_central[:, pnt], label='central')
    plt.plot(space_time_array_backward[:, pnt], label='backward')
    plt.legend(loc='best')
    plt.xlabel("time")
    plt.ylabel("concentration")
    plt.grid()
    plt.show()

    """
    #central scheme
    space_time_array_central_forward = space_time(0.5, 0, v, D, c0, col_dx, col_nx, time_dt, k, n_timesteps) 
    space_time_array_central_central = space_time(0.5, 0.5, v, D, c0, col_dx, col_nx, time_dt, k, n_timesteps) 
    space_time_array_central_backward = space_time(0.5, 1, v, D, c0, col_dx, col_nx, time_dt, k, n_timesteps) 
    #print("space_time_array_central_central",space_time_array_central_central)

    #explicit scheme
    space_time_array_explicit_forward = space_time(0, 0, v, D, c0, col_dx, col_nx, time_dt, k, n_timesteps) 
    space_time_array_explicit_central = space_time(0, 0.5, v, D, c0, col_dx, col_nx, time_dt, k, n_timesteps) 
    space_time_array_explicit_backward = space_time(0, 1, v, D, c0, col_dx, col_nx, time_dt, k, n_timesteps) 

    print("sapce time implicit", space_time_array_implicit_backward[time,pnt])
    print("sapce time central", space_time_array_central_backward[time,pnt])
    print("sapce time explicit", space_time_array_explicit_backward[time,pnt])
    print("sapce time analytical", space_time_array_analytical[time,pnt])


    #Visualization

    if th == 1: #implicit scheme

        plt.subplot(211)
        plt.title("c/c0 vs. distance (implicit)")
        #plt.plot(space_time_array_implicit_forward[time, :],label='forward')
        plt.plot(space_time_array_implicit_central[time, :],label='central')
        plt.plot(space_time_array_implicit_backward[time, :],label='backward')
        plt.plot(space_time_array_analytical[time, :],label='analytical')
        plt.legend(loc='best')
        plt.xlabel("distance")
        plt.ylabel("concentration")
        plt.grid()
        #plt.show()

        plt.subplot(212)
        #plt.title("c/c0 vs. time (implicit)")
        #plt.plot(space_time_array_implicit_forward[:, pnt],label='forward')
        plt.plot(space_time_array_implicit_central[:, pnt],label='central')
        plt.plot(space_time_array_implicit_backward[:, pnt],label='backward')
        plt.plot(space_time_array_analytical[:, pnt],label='analytical')
        plt.legend(loc='best')
        plt.xlabel("time")
        plt.ylabel("concentration")
        plt.grid()
        plt.show()

    elif th == 0.5: #central scheme

        plt.subplot(211)
        plt.title("c/c0 vs. distance (Crank-Nicolson)" )
        plt.plot(space_time_array_central_forward[time, :],label='forward')
        plt.plot(space_time_array_central_central[time, :],label='central')
        plt.plot(space_time_array_central_backward[time, :],label='backward')
        plt.plot(space_time_array_analytical[time, :],label='analytical')
        plt.legend(loc='best')
        plt.xlabel("distance")
        plt.ylabel("concentration")
        plt.grid()
        #plt.show()

        plt.subplot(212)
        #plt.title("c/c0 vs. time (Crank-Nicolson)" )
        plt.plot(space_time_array_central_forward[:, pnt],label='forward')
        plt.plot(space_time_array_central_central[:, pnt],label='central')
        plt.plot(space_time_array_central_backward[:, pnt],label='backward')
        plt.plot(space_time_array_analytical[:, pnt],label='analytical')
        plt.legend(loc='best')
        plt.xlabel("time")
        plt.ylabel("concentration")
        plt.grid()
        plt.show()

    elif th == 0:

        plt.subplot(211)
        plt.title("c/c0 vs. distance (explicit)" )
        #plt.plot(space_time_array_explicit_forward[time, :],label='forward')
        plt.plot(space_time_array_explicit_central[time, :],label='central')
        plt.plot(space_time_array_explicit_backward[time, :],label='backward')
        plt.plot(space_time_array_analytical[time, :],label='analytical')
        plt.legend(loc='best')
        plt.xlabel("distance")
        plt.ylabel("concentration")
        plt.grid()
        #plt.show()

        plt.subplot(212)
        #plt.title("c/c0 vs. time (explicit)" )
        #plt.plot(space_time_array_explicit_forward[:, pnt],label='forward')
        plt.plot(space_time_array_explicit_central[:, pnt],label='central')
        plt.plot(space_time_array_explicit_backward[:, pnt],label='backward')
        plt.plot(space_time_array_analytical[:, pnt],label='analytical')
        plt.legend(loc='best')
        plt.xlabel("time")
        plt.ylabel("concentration")
        plt.grid()
        plt.show()


    print("Done! Yay!")

    return


if __name__ == '__main__':
    main()