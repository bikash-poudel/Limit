import numpy as np
import matplotlib.pyplot as plt
import os
import time

import iso_cmf as iso
import iso_cmf_Barnes_Alison as Ba

start = time.time()

dt, sim = 12, 250
p_iso, delta, X, Y = iso.run_testcases(test_cases=[1, 2, 3, 4, 5, 6],
                                       solutes=["2H", "18O"],
                                       dt=dt,  # hours
                                       sim_period=sim,  # days
                                       )
end = time.time()
print(f'time_elapsed: {end - start}')

"""

dt, sim = 12, 250
p_iso, delta, X, Y = Ba.run_testcases(solutes=['2H', "18O"],
                                      dt=dt,  # hours
                                      sim_period=sim,  # days
                                      BA=True
                                      )
"""

er, ev, ql0, qv0 = Y
theta, m_pot, ql, qv, T = X

#delta_2H, delta_18O = delta[3]['2H'], delta[3]['18O']
delta_2H, delta_18O = delta[1]['2H'], delta[1]['18O']

##################################################################
####################### Plot information #########################
##################################################################

dz = 10 # plot layers from top
days = [50, 100, 150, 200, 250]  # profile plot days

dt_step = dt / 24  # days
f = 86400 / p_iso.get_cells()[0].area  # factor ms-1 to mm day-1

depth = [l.center.z for l in p_iso.get_cells()[0].layers][:dz]
t_step = [int(d / dt_step) - 1 for d in days]
time = np.arange(dt_step, sim + dt_step, dt_step)

output_dir = r"D:\Isotope transport\Scripts\output\iso_cmf_final" \
             r"\area=1"
os.makedirs(output_dir, exist_ok=True)

########################## Iso_delta #############################
fig1 = iso.visualize(p_iso, delta, dz, Isotopologue='2H')
fig1.savefig(os.path.join(output_dir, "2H.png"), dpi=300)

fig2 = iso.visualize(p_iso, delta, dz, Isotopologue='18O')
fig2.savefig(os.path.join(output_dir, "18_O.png"), dpi=300)
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
######################### T ######################################

T_dt = [T[t] for t in t_step]
for th, l, in zip(T_dt, days):
    plt.plot(th[:dz], depth, label=str(l) + ' days')

plt.title('Temperature')
plt.xlabel('[K]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()

plt.savefig(os.path.join(output_dir, "Temp_profile.png"), dpi=300)
plt.show()
plt.close()

##################################################################
####################### TEMPORAL #################################
##################################################################

################### surface fluxes ###############################
plt.plot(time[:250], np.array(ev[:250]) * f, label='evaporatiom')
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
d1 = [d[0] for d in delta_18O][:250]
d2 = [d[1] for d in delta_18O][:250]
d3 = [d[2] for d in delta_18O][:250]
d4 = [d[3] for d in delta_18O][:250]
d5 = [d[4] for d in delta_18O][:250]

plt.plot(time[:250], np.array(d1), label='delta l1')
plt.plot(time[:250], np.array(d2), label='delta l2')
plt.plot(time[:250], np.array(d3), label='delta l3')
plt.plot(time[:250], np.array(d4), label='delta l4')
plt.plot(time[:250], np.array(d5), label='delta l5')
plt.title('delta')
plt.xlabel('[days]')
plt.ylabel('delta')
plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "delta_temporal.png"), dpi=300)
plt.show()
plt.close()
######################################################################
################### temporal Temperature ##############################
d1 = [d[0] - 273.17 for d in T][:250]
d2 = [d[1] - 273.17 for d in T][:250]
d3 = [d[2] - 273.17 for d in T][:250]
d4 = [d[3] - 273.17 for d in T][:250]
d5 = [d[4] - 273.17 for d in T][:250]

plt.plot(time[:250], np.array(d1), label='delta l1')
plt.plot(time[:250], np.array(d2), label='delta l2')
plt.plot(time[:250], np.array(d3), label='delta l3')
plt.plot(time[:250], np.array(d4), label='delta l4')
plt.plot(time[:250], np.array(d5), label='delta l5')
plt.title('Temperature')
plt.xlabel('days')
plt.ylabel('K')
plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "Temperature.png"), dpi=300)
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
# plt.plot(time, np.array(er["2H"]), label='error')

plt.title('error')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
#plt.yscale('log')
plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "error.png"), dpi=300)
plt.show()
plt.close()
""""
#################################################################
####################### Error ###################################
plt.plot(time[:250], np.array(q_ev[:250]), label='q_ev')

plt.xlabel('[days]')
plt.ylabel('[del_S]')

plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "q_ev.png"), dpi=300)
plt.show()
plt.close()

#################################################################
####################### Error ###################################
plt.plot(time[:250], np.array(qi_ev[:250]), label='qi_ev')

plt.xlabel('[days]')
plt.ylabel('[i_storage]')

plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "qi_ev.png"), dpi=300)
plt.show()
plt.close()

#################################################################
#################################################################
plt.plot(time[:], np.array(dS[:]), label='dS')

plt.xlabel('[days]')
plt.ylabel('[dS]')

plt.grid()
plt.legend()

plt.savefig(os.path.join(output_dir, "dS.png"), dpi=300)
plt.show()
plt.close()

#################################################################
"""

#################################################################
################# Plausibility tests ############################
# Example data
x = depth

c = p_iso.get_cells()[0]

q_total = c.liquid_fluxes + c.vapor_fluxes

# Example: multiple lines per subplot
# Each subplot has a list of (y-values, label, color, linestyle)
plot_data = [

    [(np.array(theta[-1][:dz]), 'theta', 'tab:blue', '-')],

    [(np.array(_ql[-1][:dz]) * f, '$q_l$', 'tab:pink', 'dotted'),
     (np.array(_qv[-1][:dz]) * f, '$q_l$', 'tab:blue', '--'),
     (np.array(q_total[:dz]) * f, '$q$', 'tab:brown', '-.')],

    [(np.array(delta[1]["2H"][-1][:dz]), 'Test 1', 'tab:blue', '-'),
     (np.array(delta[2]["2H"][-1][:dz]), 'Test 2', 'tab:green', '--'),
     (np.array(delta[3]["2H"][-1][:dz]), 'Test 3', 'tab:red', '-.'),
     (np.array(delta[4]["2H"][-1][:dz]), 'Test 4', 'tab:purple', '--'),
     (np.array(delta[5]["2H"][-1][:dz]), 'Test 5', 'tab:orange', '-.'),
     (np.array(delta[6]["2H"][-1][:dz]), 'Test 6', 'tab:brown', '-')],

    [(np.array(delta[1]["18O"][-1][:dz]), 'Test 1', 'tab:blue', '-'),
     (np.array(delta[2]["18O"][-1][:dz]), 'Test 2', 'tab:green', '--'),
     (np.array(delta[3]["18O"][-1][:dz]), 'Test 3', 'tab:red', '-.'),
     (np.array(delta[4]["18O"][-1][:dz]), 'Test 4', 'tab:purple', '--'),
     (np.array(delta[5]["18O"][-1][:dz]), 'Test 5', 'tab:orange', '-.'),
     (np.array(delta[6]["18O"][-1][:dz]), 'Test 6', 'tab:brown', '-')],
]

x_labels = ['Volumetric water content [m$^{3}$ m$^{-3}$]',
            'moisture flux [mm day${-1}$]',
            '$\delta^{2}$H (‰)',
            '$\delta^{18}$O (‰)'
            ]

# Figure setup
fig, axs = plt.subplots(2, 2, figsize=(10, 8), dpi=300)
axs = axs.flatten()

# Iterate over subplots
for i, ax in enumerate(axs):
    for y, label, color, style in plot_data[i]:
        ax.plot(y, x, color=color, linestyle=style, linewidth=1.6, label=label)

    # Axis labels and formatting
    ax.set_xlabel(x_labels[i], fontsize=10)
    ax.set_ylabel('depth [m]', fontsize=10)
    ax.tick_params(axis='both', labelsize=9)

    # Legend (individual per subplot)
    if i != 0:
        ax.legend(fontsize=8, frameon=False, loc='best')

    # Clean journal style
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.grid(False)

# Layout adjustment for publication
fig.tight_layout(pad=0.5)

plt.show()
# Save figure
fig.savefig(os.path.join(output_dir, "palusibility.png"), dpi=300, bbox_inches='tight')
plt.show()
plt.close()


#################################################################
################# Analytical Unsaturated ############################

d_2H_analytical = np.array([7.41089387e+01, 1.98308705e+01, 5.30655857e+00, 1.41998627e+00, 3.79975268e-01, 1.01677887e-01, 2.72080676e-02, 7.28062869e-03,  1.94822928e-03, 5.21328240e-04, 1.39502643e-04, 3.73296243e-05,  9.98906418e-06, 2.67298171e-06, 7.15265325e-07, 1.91398423e-07, 5.12164579e-08, 1.37050531e-08, 3.66734619e-09, 9.81348114e-10,       2.62599731e-10, 7.02692730e-11, 1.88034112e-11, 5.03161993e-12,       1.34641522e-12, 3.60288331e-13, 9.64098441e-14, 2.57983877e-14,      6.90341130e-15, 1.84728938e-15, 4.94317650e-16, 1.32274857e-16,       3.53955354e-17, 9.47151976e-18, 2.53449159e-18, 6.78206641e-19,       1.81481860e-19, 4.85628768e-20, 1.29949792e-20, 3.47733694e-21,       9.30503387e-22, 2.48994149e-22, 6.66285445e-23, 1.78291858e-23,       4.77092616e-24, 1.27665597e-24, 3.41621396e-25, 9.14147439e-26,       2.44617448e-26, 6.54573796e-27, 1.75157928e-27, 4.68706508e-28,       1.25421551e-28, 3.35616537e-29, 8.98078989e-30, 2.40317678e-30,       6.43068007e-31, 1.72079085e-31, 4.60467807e-32, 1.23216951e-32,       3.29717228e-33, 8.82292982e-34, 2.36093488e-34, 6.31764463e-35,       1.69054360e-35, 4.52373922e-36, 1.21051101e-36, 3.23921615e-37,       8.66784453e-38, 2.31943549e-38, 6.20659606e-39, 1.66082803e-39,       4.44422307e-40, 1.18923322e-40, 3.18227874e-41, 8.51548527e-42,       2.27866555e-42, 6.09749946e-43, 1.63163478e-43, 4.36610462e-44,       1.16832945e-44, 3.12634216e-45, 8.36580410e-46, 2.23861224e-46,       5.99032050e-47, 1.60295468e-47, 4.28935931e-48, 1.14779310e-48,       3.07138879e-49, 8.21875396e-50, 2.19926298e-50, 5.88502548e-51,       1.57477870e-51, 4.21396298e-52, 1.12761774e-52, 3.01740138e-53,       8.07428860e-54, 2.16060538e-54, 5.78158129e-55, 1.54709799e-55])
d_18O_analytical = np.array([29.152451627793674, 7.670062601225156, 2.018007303736505, 0.530941361193505, 0.13969162971018653, 0.03675311971782126, 0.009669812083908984, 0.0025441449992821288, 0.0006693691378080751, 0.00017611222739912692, 4.6335444656211085e-05, 1.2190939057417452e-05, 3.207458052994949e-06, 8.4388799855928e-07, 2.2202845441655067e-07, 5.841608679678271e-08, 1.5369377792664077e-08, 4.0437110167167984e-09, 1.0639076615399229e-09, 2.799160245635861e-10, 7.36464109056914e-11, 1.9376503534394283e-11, 5.097993026424311e-12, 1.3412911597460422e-12, 3.5289612321708557e-13, 9.284760648480582e-14, 2.4428372721607505e-14, 6.427148920887184e-15, 1.6909944727805557e-15, 4.449036955844579e-16, 1.1705496471507098e-16, 3.0797372331211464e-17, 8.102844204993882e-18, 2.1318729242320323e-18, 5.608996113084057e-19, 1.4757369934666838e-19, 3.882690645489962e-20, 1.02154291146158e-20, 2.687697824109586e-21, 7.071381449251128e-22, 1.860493212899745e-22, 4.894991198095548e-23, 1.2878810125884551e-23, 3.388438171719647e-24, 8.915041942027391e-25, 2.3455636136861444e-25, 6.171220171059811e-26, 1.6236591571202393e-26, 4.2718765259150815e-27, 1.1239384185182709e-27, 2.9571022499316524e-28, 7.780189352436986e-29, 2.0469818505996173e-29, 5.385646167302803e-30, 1.4169732199084696e-30, 3.7280820974232576e-31, 9.8086512362074e-32, 2.580673830655417e-32, 6.789799392240195e-33, 1.7864084658523815e-33, 4.7000728924572116e-34, 1.2365976548297746e-34, 3.2535107325347066e-35, 8.56004541604574e-36, 2.2521633874458128e-36, 5.925482491300164e-37, 1.5590051303748751e-37, 4.101770615479254e-38, 1.079183246687807e-38, 2.8393505856630137e-39, 7.470382600034098e-40, 1.965471135289934e-40, 5.171190005235166e-41, 1.3605494168857256e-41, 3.5796300540379797e-42, 9.4180712326513e-43, 2.4779115272885918e-43, 6.519429918710846e-44, 1.7152737697414963e-44, 4.512916224038478e-45, 1.187356397822085e-45, 3.12395609725601e-46, 8.219184834042711e-47, 2.1624823535611138e-47, 5.689530073705743e-48, 1.496925623753408e-48, 3.938438313922199e-49, 1.0362102235699019e-49, 2.7262877867991326e-50, 7.172912337077239e-51, 1.8872061726030233e-51, 4.965273476856867e-52, 1.3063724068882786e-52, 3.4370893636243325e-53, 9.043044105377999e-54, 2.3792412137221874e-54, 6.2598265441476445e-55, 1.64697165368582e-55, 4.333212124832043e-56, 1.1400759252152452e-56])

# d_2H_model = [72.02477348680402, 12.971688275806414, 2.3368689473626247, 0.42055405412577684, 0.07572225907348518, 0.01396935012776268, 0.0016741040909096228, 0.0009380323064966944, 0.00034619857336437576, -0.0005544120182943146, 0.0002525444047218883, -0.00010399432681662546, 0.0003183431787068258, -0.0003376453220349518, -7.90017462559689e-06, -0.0007891789488567724, 0.0004157853616959528, 0.00014768707190526698, 7.958778036609715e-05, -0.0004977714941345468, 0.000892048164002901, -9.351037499882864e-06, -0.00046528721842076237, 8.610356871940894e-05, -0.0004316747459487402, 0.0008138894160847343, -0.00041566173147877805, 0.00022162930179980833, -8.67428979756113e-05, 0.00018456945460876284, -0.0006318752719947085, 0.0005280207113766977, -4.144882193024557e-05, 1.184483755523047e-05, 2.509769281289209e-06, 0.00036514828094524887, -0.0007722900816631295, 0.0008444957539932574, -0.0003616695843922102, 0.0005373078519887997, -0.0009248876204903667, 3.3117596887066725e-05, 0.00038235875998715585, -0.0003104047117119535, 0.0006906236673653865, -0.0007896815462649087, 0.00013547753963116804, 0.00031179877546527734, -1.9840927123482288e-05, -0.00039994456191294603]
d_2H_model = np.array(delta[6]["2H"][-1][:dz])
d_18O_model = np.array(delta[6]["18O"][-1][:dz])

# ---------- Unified publication style ----------
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.linewidth': 1.0,
    'lines.linewidth': 1.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'figure.figsize': (6.8, 4.0),  # slightly wider for visual balance
})

# ---------- Create subplots (no shared y-axis) ----------
fig, axes = plt.subplots(1, 2)

# Common colors and styles
col_analytical = 'tab:blue'
col_model = 'tab:red'
ls_analytical = '--'
ls_model = '-'

# ---------- δ²H ----------
axes[0].plot(depth, d_2H_analytical[:len(depth)], color=col_analytical, linestyle=ls_analytical, label='Analytical')
axes[0].plot(depth, d_2H_model, color=col_model, linestyle=ls_model, label='Numerical')

axes[0].set_xlabel(r'$\delta^{2}$H (‰)')
axes[0].set_ylabel('Depth (m)')
axes[0].legend(frameon=False, loc='best')

# ---------- δ¹⁸O ----------
axes[1].plot(depth, d_18O_analytical[:len(depth)],  color=col_analytical, linestyle=ls_analytical, label='Analytical')
axes[1].plot( depth, d_18O_model, color=col_model, linestyle=ls_model, label='Numerical')

axes[1].set_xlabel(r'$\delta^{18}$O (‰)')
axes[1].set_ylabel('Depth (m)')
axes[1].legend(frameon=False, loc='best')

# 🔹 Turn off top/right spines and grid for all subplots
for ax in axes.flat:
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.grid(False)

# ---------- Layout control ----------
fig.tight_layout(pad=0.8, w_pad=4.0)  # w_pad increased for spacing
plt.savefig(os.path.join(output_dir, "analytical_saturated.png"), bbox_inches='tight')
plt.show()

#################################################################
################# Spatio temporal roboustness analysis ##########

# 1: 10 layers, top layer = 2 cm
# 2 : 20 layers, top layer = 0.8 cm
# 3 : 100 layers, top layer = 1 cm

# d_2H_model = [72.02477348680402, 12.971688275806414, 2.3368689473626247, 0.42055405412577684, 0.07572225907348518, 0.01396935012776268, 0.0016741040909096228, 0.0009380323064966944, 0.00034619857336437576, -0.0005544120182943146, 0.0002525444047218883, -0.00010399432681662546, 0.0003183431787068258, -0.0003376453220349518, -7.90017462559689e-06, -0.0007891789488567724, 0.0004157853616959528, 0.00014768707190526698, 7.958778036609715e-05, -0.0004977714941345468, 0.000892048164002901, -9.351037499882864e-06, -0.00046528721842076237, 8.610356871940894e-05, -0.0004316747459487402, 0.0008138894160847343, -0.00041566173147877805, 0.00022162930179980833, -8.67428979756113e-05, 0.00018456945460876284, -0.0006318752719947085, 0.0005280207113766977, -4.144882193024557e-05, 1.184483755523047e-05, 2.509769281289209e-06, 0.00036514828094524887, -0.0007722900816631295, 0.0008444957539932574, -0.0003616695843922102, 0.0005373078519887997, -0.0009248876204903667, 3.3117596887066725e-05, 0.00038235875998715585, -0.0003104047117119535, 0.0006906236673653865, -0.0007896815462649087, 0.00013547753963116804, 0.00031179877546527734, -1.9840927123482288e-05, -0.00039994456191294603]
d_2H_model_1 = [-9.797188538656187, 40.48519023945474, 18.532047779164174, 7.283689101671209, 3.0383713899293507, 1.339876031465348, 0.632253232817126, 0.22894070559686952, 0.03901488637070294, -0.05825375897205998]
depth_1 = [-0.01, -0.035, -0.07500000000000001, -0.125, -0.175, -0.225, -0.275, -0.35, -0.5, -0.8]

d_2H_model_2 = [-5.288368096356355, 35.730868464182166, 26.232135798849534, 17.351865256761243, 11.438344085965113, 8.297447771732136, 6.055567920304261, 4.439755068431417, 2.7961182631544013, 0.8778205185844268, 0.0916161261090842, -0.049652813669753826, -0.07656364963259943, -0.07965552597877235, -0.07636580504299939, -0.07054126247352244, -0.06005644586026815, -0.025816376977294553, 0.006687740372379736, 3.2128629978700474e-05]
depth_2 = [-0.004, -0.012, -0.0205, -0.0325, -0.045, -0.055, -0.065, -0.07500000000000001, -0.09, -0.125, -0.175, -0.225, -0.275, -0.32499999999999996, -0.375, -0.42500000000000004, -0.525, -0.7, -0.875, -0.975]

d_2H_model = np.array(delta[6]["2H"][-1][:dz])

# ---------- Unified publication style ----------
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.linewidth': 1.0,
    'lines.linewidth': 1.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'figure.figsize': (6.8, 5.0),  # slightly wider for visual balance
})

# ---------- Create subplots (no shared y-axis) ----------
fig, axes = plt.subplots()

# Common colors and styles
col_analytical = 'tab:blue'
col_model = 'tab:red'
ls_analytical = '--'
ls_model = '-'

# ---------- δ²H ----------
axes.plot(d_2H_model[:50], depth[:50], linestyle=ls_model, label='100 layers, top layer = 1 cm')
axes.plot(d_2H_model_2, depth_2, linestyle=ls_model, label='20 layers, top layer = 2 cm')
axes.plot(d_2H_model_1, depth_1, linestyle=ls_model, label='10 layers, top layer = 0.8 cm')

axes.set_xlabel(r'$\delta^{2}$H (‰)')
axes.set_ylabel('Depth (m)')
axes.legend(frameon=False, loc='best')

# 🔹 Turn off top/right spines and grid for all subplots

for spine in ['top', 'right']:
    axes.spines[spine].set_visible(False)
axes.grid(False)

# ---------- Layout control ----------
fig.tight_layout(pad=0.8, w_pad=4.0)  # w_pad increased for spacing
plt.savefig(os.path.join(output_dir, "spatio_temporal-analysis.png"), bbox_inches='tight')
plt.show()







