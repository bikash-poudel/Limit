'''
Created on 03.06.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

import os
import numpy as np

from _sli_ import Sli, iso_project, iso_storages, iso_delta


class hydrus():

    def __init__(self, path):
        self.__path = path
        self.__scaler = self.__read_sli(self.__path + '/scaler.txt')
        self.__in_soil = self.__read_sli(self.__path + '/in_soil.txt')
        self.__in_variables = self.__read_sli(self.__path + '/in_variables.txt')
        self.__in_fluxes = self.__read_sli(self.__path + '/in_fluxes.txt')
        self.__in_parameter = self.__read_sli(self.__path + '/in_parameter.txt')
        self.__in_iso = self.__read_sli(self.__path + '/in_iso.txt')
        self.__inout_water = self.__read_sli(self.__path + '/inout_water.txt')
        self.__inout_iso = self.__read_sli(self.__path + '/inout_iso.txt')
        self.__out_iso = self.__read_sli(self.__path + '/out_iso.txt')


