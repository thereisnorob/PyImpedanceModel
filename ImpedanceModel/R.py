from __future__ import annotations
import numpy as np
import numpy.typing as npt
from .ImpedanceModelElement import ImpedanceModelElement

class R(ImpedanceModelElement):
    """
    Equivalent circuit element: Ideal Ohmic resistance
    
    @version:   AB-20230311
                (AA-20221224)
    @author:    Robert Leonhardt <mail@robertleonhardt.de>
    """
    
    # Default min/max values
    _min_value_list: List[float] = [1e-4]
    _max_value_list: List[float] = [1e4]
    
    
    def __init__(self, R_Ohm: float = 0.01):
        """
        Equivalent circuit element: Ideal Ohmic resistance
        
        Model function for Ohmic resistances (electrolyte, electrodes etc.)
        Appearance: Single point at R + 0i
        Model:      Z = R

        Args:
            R_Ohm (float): Value of Ohmic resistance in Ohm
        """
        self.set_parameters(R_Ohm)
        
    
    def evaluate(self, frequency_Hz: npt.ArrayLike) -> npt.ArrayLike:
        return np.zeros(shape = len(frequency_Hz)) + self.R_Ohm + 0j
    
    
    @property 
    def R_Ohm(self):
        return self._parameter_list[0]
    
    @R_Ohm.setter
    def R_Ohm(self, value):
        self._parameter_list[0] = value
        