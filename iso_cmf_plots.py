import numpy as np
import matplotlib.pyplot as plt
import os

import iso_cmf as iso
import iso_cmf_Barnes_Alison as Ba


dt, sim = 1, 250
p_iso, delta, X, Y = iso.run_testcases(test_cases=[1, 2, 3, 4, 5, 6],
                                       solutes=["18O"],
                                       dt=dt,               # hours
                                       sim_period=sim,      # days
                                       )
"""

dt, sim = 1, 250
p_iso, delta, X, Y = Ba.run_testcases(solutes=['2H', "18O"],
                                      dt=dt,  # hours
                                      sim_period=sim,  # days
                                      BA=True
                                      )
"""

theta, m_pot, ql, qv = X
er, ev, ql0, qv0 = Y
# delta_2H, delta_18O = delta[1]['2H'], delta[1]['18O']
delta_2H, delta_18O = delta[6]['2H'], delta[6]['18O']

##################################################################
####################### Plot information #########################
##################################################################

dz = 10  # plot layers from top
days = [50, 100, 150, 200, 250]  # profile plot days

dt_step = dt / 24  # days
f = 86400 * 1000  # factor ms-1 to mm day-1

depth = [-l.center.z for l in p_iso.get_cells()[0].layers][:dz]
t_step = [int(d / dt_step) - 1 for d in days]
time = np.arange(dt_step, sim + dt_step, dt_step)

output_dir = r"D:\Isotope transport\Scripts\output" \
             r"\iso_cmf_final\Barnes_Alison_unsaturated"
os.makedirs(output_dir, exist_ok=True)

########################## Iso_delta #############################
fig1 = iso.visualize(p_iso, delta, dz, Isotopologue='2H')
fig1.savefig(os.path.join(output_dir, "2H.png"), dpi=300)

fig2 = iso.visualize(p_iso, delta, dz, Isotopologue='18O')
fig2.savefig(os.path.join(output_dir, "18O.png"), dpi=300)
##################################################################

##################################################################
####################### PROFILES #################################
##################################################################

######################### delta profile ##########################
_delta = [delta_18O[t] for t in t_step]
for d, l in zip(_delta, days):
    plt.plot(d[:dz], depth, label=str(l) + ' days')

plt.title('delta')
plt.xlabel('delta')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()

plt.savefig(os.path.join(output_dir, "delta.png"), dpi=300)
plt.show()
plt.close()
##################################################################
##################### theta ######################################
theta_dt = [theta[t] for t in t_step]
for th, l, in zip(theta_dt, days):
    plt.plot(th[:dz], depth, label=str(l) + ' days')

plt.title('theta')
plt.xlabel('[m3 / m3 ]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()

plt.savefig(os.path.join(output_dir, "theta_profile.png"), dpi=300)
plt.show()
plt.close()
##################################################################
##################### matrix potential ###########################
pot = [m_pot[t] for t in t_step]
for pt, l, in zip(pot, days):
    plt.plot(pt[:dz], depth, label=str(l) + ' days')

plt.title('matrix potential')
plt.xlabel('[m]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()

plt.savefig(os.path.join(output_dir, "matrix pot.png"), dpi=300)
plt.show()
plt.close()
##################################################################
####################### ql profile ###############################
_ql = [ql[t] for t in t_step]
_qv = [qv[t] for t in t_step]
for q_l, l in zip(_ql, days):
    plt.plot(np.array(q_l[:dz]) * f, depth, label=str(l) + ' days')

plt.title('ql')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()

plt.savefig(os.path.join(output_dir, "ql_profile.png"), dpi=300)
plt.show()
plt.close()
##################################################################
##################### qv profile #################################
for q_v, l in zip(_qv, days):
    plt.plot(np.array(q_v[:dz]) * f, depth, label=str(l) + ' days')

plt.title('qv')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()

plt.savefig(os.path.join(output_dir, "qv_profile.png"), dpi=300)
plt.show()
plt.close()
##################################################################
######################### ql qv ##################################
plt.plot(np.array(_ql[-1][:dz]) * f, depth, label='ql')
plt.plot(np.array(_qv[-1][:dz]) * f, depth, label='qv')

plt.title('ql_qv after 250 days')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()

plt.savefig(os.path.join(output_dir, "ql_qv_profile.png"), dpi=300)
plt.show()
plt.close()
##################################################################

##################################################################
####################### TEMPORAL #################################
##################################################################

################### surface fluxes ###############################
plt.plot(time, np.array(ev) * f, label='evaporatiom')
# plt.plot(time, et, label='et_pot_sli')
# plt.plot(time, np.array(ql0) * f, label='ql0')
# plt.plot(time, np.array(qv0) * f, label='qv0')

plt.title('qevap')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.legend()
plt.grid()

plt.savefig(os.path.join(output_dir, "surface_fluxes.png"), dpi=300)
plt.show()
plt.close()
##################################################################
################### temporatl theta ##############################
th1 = [t[0] for t in theta]
th2 = [t[1] for t in theta]
th3 = [t[2] for t in theta]
th4 = [t[3] for t in theta]
th5 = [t[4] for t in theta]

plt.plot(time, th1, label='theta l1')
plt.plot(time, th2, label='theta l2')
plt.plot(time, th3, label='theta l3')
plt.plot(time, th4, label='theta l4')
plt.plot(time, th5, label='theta l5')

plt.title('theta')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "theta_temporal.png"), dpi=300)
plt.show()
plt.close()
##################################################################
################### temporal ql qv ###############################
ql1 = [q[0] * f for q in ql]
ql2 = [q[1] * f for q in ql]
ql3 = [q[2] * f for q in ql]
ql4 = [q[3] * f for q in ql]

qv1 = [q[0] * f for q in qv]
qv2 = [q[1] * f for q in qv]
qv3 = [q[2] * f for q in qv]
qv4 = [q[3] * f for q in qv]

plt.plot(time, ql1, label='ql_l1')
plt.plot(time, ql2, label='ql_l2')
plt.plot(time, ql3, label='ql_l3')
plt.plot(time, ql4, label='ql_l4')
plt.plot(time, qv1, linestyle='--', label='qv_l1')
plt.plot(time, qv2, linestyle='--', label='qv_l2')
plt.plot(time, qv3, linestyle='--', label='qv_l3')
plt.plot(time, qv4, linestyle='--', label='qv_l4')

plt.title('ql_qv')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "ql_qv_temporal.png"), dpi=300)
plt.show()
plt.close()
#################################################################
################### temporal delta ##############################
d1 = [d[0] for d in delta_18O]
d2 = [d[1] for d in delta_18O]
d3 = [d[2] for d in delta_18O]
d4 = [d[3] for d in delta_18O]
d5 = [d[4] for d in delta_18O]

plt.plot(time, np.array(d1), label='delta l1')
plt.plot(time, np.array(d2), label='delta l2')
plt.plot(time, np.array(d3), label='delta l3')
plt.plot(time, np.array(d4), label='delta l4')
plt.plot(time, np.array(d5), label='delta l5')
plt.title('delta')
plt.xlabel('[days]')
plt.ylabel('delta')
plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "delta_temporal.png"), dpi=300)
plt.show()
plt.close()
##################################################################
############ cumulative evaporation ##############################
"""
plt.plot(time, np.array(ev) * f, label='evaporatiom')
plt.plot(time, np.cumsum(ev) * f, label='cum evaporatiom')

plt.title('cum Evaporation')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.yscale('log')
plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "cum_ev.png"), dpi=300)
plt.show()
plt.close()
"""
#################################################################
####################### Error ###################################
plt.plot(time, np.array(er["18O"]), label='error')

plt.title('error')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.yscale('log')
plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "error.png"), dpi=300)
plt.show()
plt.close()
#################################################################
#################################################################




