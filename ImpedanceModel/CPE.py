from __future__ import annotations
import numpy as np
import numpy.typing as npt
from .ImpedanceModelElement import ImpedanceModelElement

class CPE(ImpedanceModelElement):
    """
    Equivalent circuit element: Non-ideal constant phase element
    
    @version:   AB-20230311
                (AA-20221224)
    @author:    Robert Leonhardt <mail@robertleonhardt.de>
    """
    
    # Default min/max values
    _min_value_list: List[float] = [1e-6, 0]
    _max_value_list: List[float] = [1e6, 1]
    
    
    def __init__(self, Q0_Ohm_p_s_n: float = 30, n_1: float = 0.95):
        """
        Equivalent circuit element: Non-ideal constant phase element
        
        Model function for constant phase elements (CPE, non-ideal capacitors)
        Appearance: Line with slope of n and length ~ 1/Q0
        Model:      Z = 1 / ((jw)^n * Q0)

        Args:
            Q0_Ohm_p_s_n (float): Capacitance value in Ohm/s^n
            n_1 (float): Phase angle (0 ... 1) whereas 0 = point (ideal resistor) and 1 = vertical line (ideal capacitor)
        """
        self.set_parameters(Q0_Ohm_p_s_n, n_1)
        
    
    def evaluate(self, frequency_Hz: npt.ArrayLike) -> npt.ArrayLike:
        return 1 / (self.Q0_Ohm_p_s_n * (1j * 2 * np.pi * frequency_Hz) ** self.n_1)
    
    
    @property 
    def Q0_Ohm_p_s_n(self):
        return self._parameter_list[0]
    
    @Q0_Ohm_p_s_n.setter
    def Q0_Ohm_p_s_n(self, value):
        self._parameter_list[0] = value
        
    @property 
    def n_1(self):
        return self._parameter_list[1]
    
    @n_1.setter
    def n_1(self, value):
        self._parameter_list[1] = value