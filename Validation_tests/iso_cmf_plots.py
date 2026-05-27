
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
                                       check_mass=False
                                       )

"""

p_iso, delta, X, Y = Ba.run_testcases(solutes=['2H', '18O'],
                                      dt=dt,  # hours
                                      sim_period=sim,  # days
                                      BA=True
                                      )
"""

end = time.time()
print(f'time_elapsed: {end - start}')

##################################################################
#----------------------- Plot information -----------------------#
##################################################################

dz = 10  # plot layers from top
# output_dir =
# os.makedirs(output_dir, exist_ok=True)

########################## Iso_delta #############################
fig1 = iso.visualize(p_iso, delta, dz, Isotopologue='2H')
# fig1.savefig(os.path.join(output_dir, "2H.png"), dpi=300)

fig2 = iso.visualize(p_iso, delta, dz, Isotopologue='18O')
# fig2.savefig(os.path.join(output_dir, "18_O.png"), dpi=300)
##################################################################


