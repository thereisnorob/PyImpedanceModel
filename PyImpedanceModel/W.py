from __future__ import annotations
import numpy as np
import numpy.typing as npt
from .ImpedanceModelElement import ImpedanceModelElement

class W(ImpedanceModelElement):
    """
    Equivalent circuit element: Ideal Warburg impedance (semi infinite)
    
    @version:   AB-20230311
                (AA-20221224)
    @author:    Robert Leonhardt <mail@robertleonhardt.de>
    """
    
    # Default min/max values
    _min_value_list: List[float] = [1e-9]
    _max_value_list: List[float] = [1e9]
    
    
    def __init__(self, Aw_Ohm_p_s_sqrt: float = 1e-6):
        """
        Equivalent circuit element: Ideal Warburg impedance (semi infinite)
        
        Model function for default Warburg diffusion
        Appearance: Line with slope of 45°
        Model:      Z = Aw / sqrt(w) + Aw / (j * sqrt(w))

        Args:
            Aw_Ohm_p_s_sqrt (float): Warburg coefficient in Ohm/s^(1/2)
        """
        self.set_parameters(Aw_Ohm_p_s_sqrt)
        
    
    def evaluate(self, frequency_Hz: npt.ArrayLike) -> npt.ArrayLike:
        return self.Aw_Ohm_p_s_sqrt / np.sqrt(2 * np.pi * frequency_Hz) + self.Aw_Ohm_p_s_sqrt / (1j * np.sqrt(2 * np.pi * frequency_Hz))
    
    
    @property 
    def Aw_Ohm_p_s_sqrt(self):
        return self._parameter_list[0]
    
    @Aw_Ohm_p_s_sqrt.setter
    def Aw_Ohm_p_s_sqrt(self, value):
        self._parameter_list[0] = value