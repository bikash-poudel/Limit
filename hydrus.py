import os
import subprocess

import numpy as np


class hydrus(object):

    def __init__(self, path):
        self.__path = path
        self.__profile = self.__read_hydrus(self.__path + '/Profile.out')
        self.__Nod_inf = self.__read_hydrus(self.__path + '/Nod_Inf.out')
        self.__Nod_Inf_V = self.__read_hydrus(self.__path + '/Nod_Inf_V.out')
        self.__T_level = self.__read_hydrus(self.__path + '/T_Level.out')
        self.__obs_node = self.__read_hydrus(self.__path + '/Obs_Node.out')
        self.__Atm = self.__read_hydrus(self.__path + '/ATMOSPH.IN')

    def get_path(self):
        return self.__path

    def __read_hydrus(self, path):
        # Read the text file
        with open(path, 'r') as file:
            # Read all lines into a list, excluding empty lines
            lines = [line.strip() for line in file if line.strip()]

        # Split each line by spaces and flatten the resulting list
        data = [line.split() for line in lines]

        return data

    def extract_nodes(self, data, time):

        data_blocks = {}
        current_time = None
        current_block = []

        for line in data:
            if "Time:" in line:
                # Store the previous block if it exists
                if current_time is not None:
                    data_blocks[current_time] = current_block

                # Extract the new time
                current_time = float(line[1])
                current_block = []
            elif current_time is not None:
                # Append the line to the current data block
                current_block.append(line)

        # Add the last block if it wasn't added
        if current_time is not None:
            data_blocks[current_time] = current_block

        # Return the data for the specified time
        return data_blocks.get(time, None)


    @property
    def time(self):

        t = []
        for line in self.__Nod_inf:
            if "Time:" in line:
                t.append(round(float(line[1]), 4))

        return t[1:]  # Ignoring the time stamp

    @property
    def dt(self):
        dt = np.diff(np.array(self.time))
        return dt


    ###### Profile #######
    @property
    def profile(self):
        return self.__profile[6:-1]

    @property
    def depth(self):
        d = list(zip(*self.profile))[1]
        return [float(x) for x in d][1:]  # covert [cm] to [m]

    @property
    def thetar(self):
        return list(map(float, list(zip(*self.profile))[2]))[1:]

    @property
    def theta_sat(self):
        return list(map(float, list(zip(*self.profile))[3]))[1:]

    @property
    def Ks(self):

        K = list(zip(*self.profile))[5]

        return [float(x) for x in K][1:]  # covert [m s-1]

    ###### Nod_inf #######
    @property
    def nod_inf(self):
        return self.__Nod_inf

    def node_output(self, time):
        return self.extract_nodes(self.nod_inf, time)[2:-1]

    def head(self, time):

        h = list(zip(*self.node_output(time)))[2]
        ht = np.array([float(x) for x in h])
        return (ht[:-1] + ht[1:]) / 2

    def theta(self, time):
        th = list(map(float, list(zip(*self.node_output(time)))[3]))
        t_h = np.array(th)
        return (t_h[:-1] + t_h[1:]) / 2

    def transpiration(self, time):  # TODO: check rtanspiratin before implementing
        tr = list(map(float, list(zip(*self.node_output(time)))[7]))
        t_r = np.array(tr)
        return (t_r[:-1] + t_r[1:]) / 2

    def T(self, time):
        T = list(map(float, list(zip(*self.node_output(time)))[10]))
        Tmp = np.array(T)
        return (Tmp[:-1] + Tmp[1:]) / 2

    ###### Nod_inf_v #######
    @property
    def node_inf_v(self):
        return self.__Nod_Inf_V

    def node_v_output(self, time):
        return self.extract_nodes(self.node_inf_v, time)[2:-1]

    def cond_iso_thermal_liquid(self, time):
        """"Isothermal hydraulic conductivity of liquid phase"""
        cond = list(zip(*self.node_v_output(time)))[2]

        return [float(x) for x in cond]

    def cond_thermal_liquid(self, time):
        """"Thermal hydraulic conductivity of liquid phase"""
        return list(map(float, list(zip(*self.node_v_output(time)))[3]))

    def cond_iso_thermal_vapor(self, time):
        """"Isothermal hydraulic conductivity of vapor phase"""

        cond = list(zip(*self.node_v_output(time)))[4]
        return [float(x) for x in cond]

    def cond_thermal_vapor(self, time):
        """"Thermal hydraulic conductivity of vapor phase"""
        return list(map(float, list(zip(*self.node_v_output(time)))[5]))

    def ql(self, time):
        q = list(zip(*self.node_v_output(time)))[6]
        return q

    def q_liquid(self, time):
        """"Returns total liquid flow cm/day"""
        return [float(x) for x in self.ql(time)]

    def ql_surface(self, time):
        """Returns top surface flow cm/day"""
        return float(self.ql(time)[0])

    def qv(self, time):
        q = list(zip(*self.node_v_output(time)))[7]
        return q

    def q_vapor(self, time):
        """"Returns total vapor flow cm/day"""
        return [float(x) for x in self.qv(time)]

    def qv_surface(self, time):
        """"Returns total vapor flow cm/day"""
        return float(self.qv(time)[0])

    def q_total(self, time):
        """"Returns total water flow cm/day"""

        q = list(zip(*self.node_v_output(time)))[8]
        return [float(x) for x in q]

    def q_vapor_isothermal(self, time):
        """"Returns isothermal vapor flow cm/day"""

        q = list(zip(*self.node_v_output(time)))[9]
        return [float(x) for x in q]

    def q_vapor_thermal(self, time):
        """"Returns thermal vapor flow cm/day"""

        q = list(zip(*self.node_v_output(time)))[10]
        return [float(x) for x in q]

    ###### T_Level ######
    @property
    def T_level(self):

        t_level = np.array([t for t in self.__T_level[7:-1] if round(float(t[0]), 2) in self.time], dtype=float)

        return t_level

    @property
    def q_evap_all(self):

        t_level = self.T_level
        time = t_level[:, 0]
        ev = t_level[:, 18]

        dt = np.diff(time)
        d_ev = np.diff(ev)
        dt = np.insert(dt, 0, time[0])
        d_ev = np.insert(d_ev, 0, ev[0])

        t_level[:, 20] = dt
        t_level[:, 21] = d_ev

        return {round(qev[0], 2): qev for qev in t_level}

    def q_evap(self, time):

        return self.q_evap_all[time][21] / self.q_evap_all[time][20]  # cm day-1 into m s-1


    ###### Atmosphere #######

    @property
    def atm(self):
        return self.__Atm[9]

    def prec(self):
        return float(self.atm[1])

    def pot_evaporation(self):
        return float(self.atm[2])

    def bot_flux(self):
        return float(self.atm[5])

    def T_s(self):
        return float(self.atm[8])

    def T_bottom(self):
        return float(self.atm[9])

    ###### Obs node ######
    @property
    def obs_node(self):
        return self.__obs_node[7:-1]

    def T_surface(self, time):

        Ts = None
        if round(time, 2) in self.time:

            for t in self.obs_node:
                if round(float(t[0]), 2) == time:
                    Ts = float(t[3])
        else:
            raise 'time parameter is not in self.time'

        return Ts


def run_hydrus(path_H1D_files, output_path):
    #write_selector_in_file(output_path)
    print('write SELECTOR.IN . done!')
    run_hydrus = subprocess.Popen([str(path_H1D_files), str(output_path)])
    # run_hydrus.communicate()


def write_selector_in_file(path):
    """
    Writes the selector file for HYDRUS-1D, including time points as fractions of days over 250 days.
    """

    # Calculate time points as fractions of days

    time_points = [hour for hour in range(1, 250 * 1 + 1)]  # 1/24 to 6000/24 days

    # Format TPrint list to match the required format (10 numbers per line)
    tprint_str = ""
    for i in range(0, len(time_points), 6):
        #tprint_str += " ".join(f"{time_point:10.05f}" for time_point in time_points[i:i + 6]) + "\n"
        tprint_str += " ".join(f"{time_points:10.00f}" for time_points in time_points[i:i + 6]) + "\n"
    file_content = f"""Pcp_File_Version=4
*** BLOCK A: BASIC INFORMATION *****************************************
Heading
Welcome to HYDRUS-1D
LUnit  TUnit  MUnit  (indicated units are obligatory for all input data)
m
sec
mmol
lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse
 t     f     t      f     f     f      f     t       t         t         f
lSnow  lHP1   lMeteo  lVapor lActiveU lFluxes lIrrig  lDummy  lDummy  lDummy
 f       f       t       t       f       f       f       f       f       f
NMat    NLay  CosAlpha
  1       1       1
*** BLOCK B: WATER FLOW INFORMATION ************************************
MaxIt   TolTh   TolH       (maximum number of iterations and tolerances)
  10    0.001   0.01
TopInf WLayer KodTop InitCond
 t     t      -1       t
BotInf qGWLF FreeD SeepF KodBot DrainF  hSeep
 f     f     f     f     -1      f      0
         rTop         rBot        rRoot
           0            0            0
    hTab1   hTabN
    1e-008     100
    Model   Hysteresis
      0          0
   thr     ths    Alfa      n         Ks       l
   0.01    0.35    5.18    2.22  1.23e-007    0.67 
*** BLOCK C: TIME INFORMATION ******************************************
        dt       dtMin       dtMax     DMul    DMul2  ItMin ItMax  MPL
       86.4       0.864      432000     1.3     0.7     3     7    {250}
      tInit        tMax
          0       {450000}
  lPrintD  nPrintSteps tPrintInterval lEnter
     f           1         86400       t
TPrint(1),TPrint(2),...,TPrint(MPL)
{tprint_str}*** BLOCK E: HEAT TRANSPORT INFORMATION *********************************************************
    Qn      Qo    Disper.    B1          B2          B3          Cn          Co           Cw
    0.6    0.01       0   -0.197001      -0.962     2.52101   1.92e+006 2.50999e+006   4.18e+006 
      tAmpl     tPeriod    Campbell   MeltConst  lDummy  lDummy  lDummy  lDummy  lDummy
          5       86400          0 4.97685e-008       f       f       f       f       f
      kTopT       TTop      kBotT       TBot
          1         20           1         20
*** END OF INPUT FILE 'SELECTOR.IN' ************************************
    """

    # Write the file content to 'SELECTOR.IN'
    f_path = path + '\SELECTOR.IN'
    with open(f_path, 'w') as writer:
        writer.write(file_content)


path_H1D_files = 'D:\Hydrus\sli_vapor'
path_H1D = 'C:\Program Files (x86)\PC-Progress\Hydrus-1D 4.xx\H1D_CALC.EXE'

# write_selector_in_file(path_H1D_files
