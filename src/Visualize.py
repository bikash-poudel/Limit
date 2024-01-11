
import numpy as np
import matplotlib.pyplot as plt


class Visualize(object):

    def __init__(self, layers, time):
        self.__layers = layers
        self.time = time

    def breakthrough(self, C, solute_i, print_steps, label_time='days'):

        time_steps = np.arange(0, self.time.time_steps) * self.time.dt * \
                     self.time.total_seconds() / self.time.total_seconds_per_day()

        # Concentration vs time
        fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
        for index, layer in enumerate(self.__layers[::print_steps]):

            value_index = index * print_steps  # index value for corresponding layer
            ax.plot(time_steps, C[value_index], label="{:.2f}m".format(layer.center))

        ax.set_xlabel('t [{}]'.format(label_time))
        ax.set_ylabel('delta')
        ax.set_title('conc breakthrough [{}]'.format(solute_i))
        ax.legend()

        # timestamp1 = datetime.now().strftime("%Y%m%d-%H%M%S")
        # fig1.savefig( 'D:\\Isotope transport\\Scripts\\output\\{}.png'.format(timestamp1))
        plt.show()
        return ax

    def profile(self, C, solute_i, print_time_steps, label_depth='m'):

        layers_center = np.array([layer.center for layer in self.__layers])

        # Concentration vs depth
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

        # timestamp2 = datetime.now().strftime("%Y%m%d-%H%M%S")
        # fig2.savefig( 'D:\\Isotope transport\\Scripts\\output\\{}.png'.format(timestamp2))
        plt.show()
        return ax

