
import iso_cmf as iso

import numpy as np
import matplotlib.pyplot as plt


p_iso, delta, X, Y = iso.run_testcases(test_cases=[1, 2, 3, 4, 5, 6])
ev, theta, ql, qv = X
T, m_pot, ql0, qv0 = Y

delta_2H_1, delta_18O_1 = delta[1]['2H'], delta[1]['18O']
delta_2H_6, delta_18O_6 = delta[1]['2H'], delta[1]['18O']

dz = 10  # plot layers from top

########################## Iso_delta #############################
iso.visualize(p_iso, delta, dz, Isotopologue='2H')
iso.visualize(p_iso, delta, dz, Isotopologue='18O')
##################################################################

f = 86400 * 1000  # factor ms-1 to mm day-1
depth = [-l.center.z for l in p_iso.get_cells()[0].layers][:dz]
days = [50, 100, 150, 200, 250]
t_step = [(d * 2)-1 for d in days]
time = np.arange(0.5, 250.5, 0.5)

_ql = [ql[t] for t in t_step]
_qv = [qv[t] for t in t_step]

############################################################
####################### PROFILES ###########################
############################################################

############ delta profile ##################
_delta_18O = [delta_18O_1[t] for t in t_step]
for d, l in zip(_delta_18O, days):
    plt.plot(d[:dz], depth, label=str(l) + ' days')

plt.title('delta')
#plt.xscale('log')
plt.xlabel('delta')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()
############################################################
##################### theta ################################
theta_dt = [theta[t] for t in t_step]
for th, l, in zip(theta_dt, days):
    plt.plot(th[:dz], depth, label=str(l) + ' days')

plt.title('theta')
plt.xlabel('[m3 / m3 ]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()
#############################################################
####################### ql profile ##########################
for q_l, l in zip(_ql, days):
    plt.plot(np.array(q_l[:dz]) * f, depth, label=str(l) + ' days')

plt.title('ql')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()
#############################################################
##################### qv profile ############################
for q_v, l in zip(_qv, days):
    plt.plot(np.array(q_v[:dz]) * f, depth, label=str(l) + ' days')

plt.title('qv')
#plt.xscale('log')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()
#############################################################
################## coupled profile ##########################
plt.plot(np.array(_ql[-1][:dz]) * f, depth, label='ql')
plt.plot(np.array(_qv[-1][:dz]) * f, depth, label='qv')
plt.title('ql_qv after 250 days')
#plt.xscale('log')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()
#############################################################

############################################################
####################### TEMPORAL ###########################
############################################################

################### surface fluxes #########################
plt.plot(time, np.array(ev) * f, label='evaporatiom')
plt.plot(time, np.array(ql0) * f, label='ql0')
plt.plot(time, np.array(qv0) * f, label='qv0')
plt.title('qevap')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()
#############################################################
################### temporatl theta #########################
th1 = [t[0] for t in theta]
th2 = [t[1] for t in theta]
th3 = [t[2] for t in theta]
th4 = [t[3] for t in theta]
th5 = [t[4] for t in theta]

# plt.plot(time, th1, label='theta l1')
plt.plot(time, th2, label='theta l2')
plt.plot(time, th3, label='theta l3')
plt.plot(time, th4, label='theta l4')
plt.plot(time, th5, label='theta l5')

plt.title('theta')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()
#############################################################
################### temporal ql qv ##########################
ql1 = [q[0] * f for q in ql]
ql2 = [q[1] * f for q in ql]
ql3 = [q[2] * f for q in ql]
ql4 = [q[3] * f for q in ql]

qv1 = [q[0] * f for q in qv]
qv2 = [q[1] * f for q in qv]
qv3 = [q[2] * f for q in qv]
qv4 = [q[3] * f for q in qv]

plt.plot(time, ql1, label='ql_l1')
plt.plot(time, ql2,  label='ql_l2')
plt.plot(time, ql3, label='ql_l3')
plt.plot(time, ql4, label='ql_l4')
plt.plot(time, qv1,  linestyle='--', label='qv_l1')
plt.plot(time, qv2,  linestyle='--', label='qv_l2')
plt.plot(time, qv3,  linestyle='--', label='qv_l3')
plt.plot(time, qv4,  linestyle='--', label='qv_l4')

plt.title('ql_qv')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()
#############################################################
################### temporal delta ##########################
d1 = [d[0] for d in delta_18O_1]
d2 = [d[1] for d in delta_18O_1]
d3 = [d[2] for d in delta_18O_1]
d4 = [d[3] for d in delta_18O_1]
d5 = [d[4] for d in delta_18O_1]

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
plt.show()
#############################################################
############ cumulative evaporation ########################
plt.plot(time, np.array(ev) * f, label='evaporatiom')
plt.plot(time, np.cumsum(ev) * f, label='cum evaporatiom')
plt.title('cum Evaporation')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()
#############################################################
#############################################################

