from __future__ import annotations
import numpy as np
import numpy.typing as npt
from .ImpedanceModelElement import ImpedanceModelElement

class FLW(ImpedanceModelElement):
    """
    Equivalent circuit element: Finite-length Warburg impedance
    
    @version:   AB-20230311
                (AA-20221224)
    @author:    Robert Leonhardt <mail@robertleonhardt.de>
    """
    
    # Default min/max values
    _min_value_list: List[float] = [1e-9, 1e-9]
    _max_value_list: List[float] = [1e9, 1e9]
    
    
    def __init__(self, Aw_Ohm_p_s_sqrt: float = 1e-7, B_1_p_s_sqrt: float = 1e-7):
        """
        Equivalent circuit element: Finite-length Warburg impedance
        
        Short (Finite-Length) Warburg element
        Appearance: Downwards-curled Warburg
        Model:      Z = Aw / sqrt(jw) * tanh(B * sqrt(jw))

        Args:
            Aw_Ohm_p_s_sqrt (float): Warburg coefficient in Ohm/s^(1/2)
            B_1_p_s_sqrt (float): B = diffusion layer thickness / diffusion coefficient in 1/sqrt(s)
        """
        self.set_parameters(Aw_Ohm_p_s_sqrt, B_1_p_s_sqrt)
        
    
    def evaluate(self, frequency_Hz: npt.ArrayLike) -> npt.ArrayLike:
        return self.Aw_Ohm_p_s_sqrt / np.sqrt(1j * 2 * np.pi * frequency_Hz) * np.tanh(self.B_1_p_s_sqrt * np.sqrt(1j * 2 * np.pi * frequency_Hz))
    
    
    @property 
    def Aw_Ohm_p_s_sqrt(self):
        return self._parameter_list[0]
    
    @Aw_Ohm_p_s_sqrt.setter
    def Aw_Ohm_p_s_sqrt(self, value):
        self._parameter_list[0] = value
        
    @property 
    def B_1_p_s_sqrt(self):
        return self._parameter_list[1]
    
    @B_1_p_s_sqrt.setter
    def B_1_p_s_sqrt(self, value):
        self._parameter_list[1] = value