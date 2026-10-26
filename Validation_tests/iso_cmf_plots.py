
import time
import importlib

import iso_cmf as iso
import Barnes_Alison_saturated_isothermal as Ba1
import Barnes_Alison_unsaturated_nonisothermal as Ba2
importlib.reload(Ba1)
importlib.reload(Ba2)


start = time.time()
dt, sim = 12, 250

run_validation_tests = "Plausibility"

if run_validation_tests == "Plausibility":

    #-------------- Plausibility Tests --------------------------------
    p_iso, delta, X, Y = iso.run_testcases(test_cases=[1, 2, 3, 4, 5, 6],
                                           solutes=["2H", "18O"],
                                           dt=dt,  # hours
                                           sim_period=sim,  # days
                                           check_mass=False
                                           )


if run_validation_tests == "BA_saturated":
    #-------------- BA saturated isthermal --------------------------------
    p_iso, delta, X, Y = Ba1.run_testcases(solutes=['2H', '18O'],
                                           bottom_potential=1.3,
                                          dt=dt,  # hours
                                          sim_period=sim,  # days
                                          BA=True
                                          )

if run_validation_tests == "BA_unsaturated":
    #-------------- BA Unsaturated non-Isothermal --------------------------
    p_iso, delta, X, Y = Ba2.run_testcases(solutes=['2H', '18O'],
                                           bottom_potential=-0.4,
                                          dt=dt,  # hours
                                          sim_period=sim,  # days
                                          BA=True
                                          )

end = time.time()
print(f'time_elapsed: {end - start}')

##################################################################
#----------------------- Plot information -----------------------#
##################################################################

dz = 50  # plot layers from top
# output_dir =
# os.makedirs(output_dir, exist_ok=True)

########################## Iso_delta #############################
fig1 = iso.visualize(p_iso, delta, dz, Isotopologue='2H')
# fig1.savefig(os.path.join(output_dir, "2H.png"), dpi=300)

fig2 = iso.visualize(p_iso, delta, dz, Isotopologue='18O')
# fig2.savefig(os.path.join(output_dir, "18_O.png"), dpi=300)
##################################################################



