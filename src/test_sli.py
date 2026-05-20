import numpy as np
import os

import matplotlib.pyplot as plt


def breakthrough(x, y, labels=None, x_label=None, y_label=None,  title=None, lgd_title=None):

    if labels is not None:
        for p, l in zip(y, labels):

            plt.plot(x, p, label=l)

        lgd = plt.legend()
        lgd.set_title(lgd_title)
    else:
        for p in y:

            plt.plot(x, p)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.show()


def profile(x, y, labels=None, x_label=None, y_label=None,  title=None, lgd_title=None):

    for p, l in zip(x, labels):

        plt.plot(p, y, label=l)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    lgd = plt.legend() if labels is not None else None
    lgd.set_title(lgd_title)
    plt.show()


def plots(plot):

    # density of liquid water
    rho_w = 1000000.  # [g/m3]
    # mole mass of water
    Mw = 18.01528  # [g/mol]

    current_dir = os.getcwd()
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir)) + '\sli_label3_von_Arndt\sli_label3'
    inpath = os.path.join(parent_dir, 'input')
    outpath = os.path.join(parent_dir, 'output')
    obspath = os.path.join(parent_dir, 'obs')

    ###############################################################################
    # load input data
    # T[degC], P[hPa], Rn[W/m2], hr[-], d18O_outlet[permil], dD_outlet[permil], prec[m/h], wetfeed[m/h], Tr[m/s], Tr_d18O[permil], Tr_dD[permil], ET[m/s], ET_d18O[permil], ET_dD[permil]
    inp = np.loadtxt(inpath + '\\bs_o_inp.csv', dtype=float, delimiter=' ')
    time = inp[:, 0].astype(float)  # jd
    rn = inp[:, 3].astype(float)  # net radiation
    # tr_inp = inp[:,8]
    # et_inp = inp[:,11]
    # d18o_tr_inp = inp[:,9]
    # d18o_et_inp = inp[:,12]
    d18o_amb = inp[:, 5]

    ###############################################################################
    # load output data
    # T[degC], P[hPa], Rn[W/m2], hr[-], d18O_outlet[permil], dD_outlet[permil], prec[m/h], wetfeed[m/h], Tr[m/s], Tr_d18O[permil], Tr_dD[permil], ET[m/s], ET_d18O[permil], ET_dD[permil]
    out = np.loadtxt(outpath + '/maren.out', dtype=float)
    step = out[:, 0].astype(int)  # timestep
    tr = out[:, 1]  # transpiration [m/s]
    ev = out[:, 2]  # evaporation [m/s]
    ev = ev * rho_w * 1.e3 / Mw  # from [m/s] to [mmol/m2s]
    Ta = out[:, 3]  # air temperature [degC]
    T0 = out[:, 4]  # surface temperature [degC]
    d18o_tr = out[:, 5]  # [permil]
    d18o_ev = out[:, 6]  # [permil]
    d18o_et = out[:, 7]  # [permil]

    out2 = np.loadtxt(outpath + '/theta.out', dtype=float)
    t = out2[1:, 0]
    h = out2[0, 1:]
    theta = out2[1:, 1:]

    out3 = np.loadtxt(outpath + '/mareniso.out', dtype=float)
    d18o_s = out3[:, 1:]

    out4 = np.loadtxt(outpath + '/soilt.out', dtype=float)
    soilt = out4[:, 1:]

    ###############################################################################
    # load obs
    obs = np.loadtxt(obspath + '/T0_o_obs.csv', delimiter=' ', skiprows=0)
    # np.place(obs, obs=='', 'nan')
    obs_jd = obs[:, 0]  # jd
    # obs_jd = date2dec(ascii=obs_ad)
    T0_obs = obs[:, 2].astype(float)  # surface temp [degC]

    ###############################################################################
    # load obs
    obs1 = np.loadtxt(obspath + '/soilt_bs_o_obs.csv', delimiter=' ', skiprows=0)
    # np.place(obs1, obs1=='', 'nan')
    obs1_jd = obs1[:, 0]  # jd
    # obs1_jd = date2dec(ascii=obs1_ad)
    soilt_obs = obs1[:, 1:5].astype(float)  # surface temp [degC]

    ###############################################################################
    # load obs
    obs2 = np.loadtxt(obspath + '/soilm_bs_o_obs.csv', delimiter=' ', skiprows=0)
    # np.place(obs2, obs2=='', 'nan')
    obs2_jd = obs2[:, 0]  # jd
    # obs2_jd = date2dec(ascii=obs2_ad)
    soilm_obs = obs2[:, 1:5].astype(float)  # surface temp [degC]

    if plot not in [0,1,2,3]:
        raise ValueError('Only numerical values: [0,1, 2,3], Please!')

    ###############################################################################
    # plot theta
    if plot == 0:
        breakthrough(t, [theta[:, 0], theta[:, 4], theta[:, 8], theta[:, 12]],
                     labels=['0', '4', '8', '12'], x_label="time", y_label='theta',
                     title='theta vs time', lgd_title='depth')

        profile([theta[0, :], theta[100, :], theta[250, :], theta[400, :], theta[500, :]], -h,
                labels=['0', '100', '250', '400', '500'], x_label="theta", y_label='depth',
                title='theta profile', lgd_title='time')

    ###############################################################################
    # plot soil isotopes with time
    elif plot == 1:
        breakthrough(t, [d18o_s[:, 0], d18o_s[:, 4], d18o_s[:, 8], d18o_s[:, 12]],
                     labels=['0', '4', '8', '12'], x_label="time", y_label='18O',
                     title=' H18O vs time', lgd_title='depth')

        profile([theta[0, :-3], d18o_s[100, :-3], d18o_s[250, :-3], d18o_s[400, :-3], d18o_s[500, :-3]], -h[0:-3],
                labels=['0', '100', '250', '400', '500'], x_label="18O", y_label='depth',
                title='18O profile', lgd_title='time')

    ###############################################################################
    # plot evaporation
    elif plot == 2:
        breakthrough(t, [ev], x_label="time", y_label='evap', title='Evaporation')

    ###############################################################################
    # plot evap isotopes
    elif plot == 3:
        breakthrough(t, [d18o_ev, d18o_amb], x_label="time", y_label='18O',
                     labels=['evap', 'amb'], title='iso_Evaporation')




