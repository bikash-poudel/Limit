
import cmf
import numpy as np
from numpy.typing import NDArray


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


