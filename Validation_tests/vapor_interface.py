
import cmf
import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt


class Vaporizer:

    """
    Adds von Neumann boundaries to every layer of a cell. The flux can
    be set as an array in m³ Water / m³ Soil / day

    !!! NOTE: If the water content becomes low and the flux is still positive
        a negative water content can occur in the layers which leads to strange
        errors and the simulation does not continue
    """

    def __init__(self, cell):
        self.cell = cell
        self.boundaries = [cmf.NeumannBoundary.create(layer) for layer in self.cell.layers]

    def __set_flux(self, vaporfluxes: NDArray[np.float64]):
        for nb, layer, vf in zip(self.boundaries, self.cell.layers, vaporfluxes):
            nb: cmf.NeumannBoundary
            layer: cmf.SoilLayer
            nb.flux = cmf.timeseries.from_scalar(-vf * layer.get_capacity())

    def __get_flux(self):
        return np.array([-nb.flux[cmf.never] / layer.get_capacity() for nb, layer in zip(self.boundaries, self.layers)])

    flux = property(__get_flux, __set_flux)


class Cmf1D(cmf.project):

    def __init__(self, layer_depth: np.ndarray):
        super().__init__()
        self.cell: cmf.Cell = self.NewCell(0, 0, 0, 1000)
        self.layers = [
            self.cell.add_layer(ld, cmf.VanGenuchtenMualem())
            for ld in layer_depth
        ]

        self.cell.install_connection(cmf.Richards)


def main():

    cmf1d = Cmf1D(np.arange(0.05, 1.05, 0.05))

    vap = Vaporizer(cmf1d.cell)
    vap.flux = [0.000001] * len(cmf1d.cell.layers)
    cmf1d.cell.layers[0].wetness = 1.0
    solver = cmf.CVodeAdams(cmf1d)

    result = []
    for t in solver.run(solver.t, solver.t + cmf.week, cmf.h):
        print(t, cmf1d.cell.layers[0].theta, cmf1d.cell.layers[-1].theta)
        result.append(cmf1d.cell.layers.wetness)


    plt.imshow(result)
    plt.show()


if __name__=='__main__':
    main()
    


