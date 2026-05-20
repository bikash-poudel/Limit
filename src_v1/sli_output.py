import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


sli_output_pth = 'D:\Isotope transport\Materials to start\Soil_Litter_Iso\sli_label3_von_Arndt\sli_label3\output'

df = pd.read_csv(sli_output_pth + '\\theta.out', delim_whitespace=True, header=None)
df = df.T
df_theta = df.iloc[1:]
df_theta.set_index(0, inplace=True)

fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
for t in df_theta.columns[1::5]:

    ax.plot(df_theta[t], - np.arange(0, 19, 1))

plt.show()

""""
fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
for column in df.columns:
    df.plot()

    plt.plot(df[column])
    plt.show()

fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
for t in range(0, self.time.time_steps, int(print_time_steps)):

    ax.plot(C[:, t], - layers_center,
            label='day {}'.format(t * self.time.dt *
                                  self.time.total_seconds() / self.time.total_seconds_per_day()))

ax.set_xlabel('delta')
ax.set_ylabel('depth [{}]'.format(label_depth))
ax.set_title('conc profile [{}]'.format(solute_i))
# plt.yscale('log')
ax.legend()

"""

