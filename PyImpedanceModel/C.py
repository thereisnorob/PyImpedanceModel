from __future__ import annotations
import numpy as np
import numpy.typing as npt
from .ImpedanceModelElement import ImpedanceModelElement

class C(ImpedanceModelElement):
    """
    Equivalent circuit element: Ideal capacitor
    
    @version:   AB-20230311
                (AA-20221224)
    @author:    Robert Leonhardt <mail@robertleonhardt.de>
    """
    
    # Default min/max values
    _min_value_list: List[float] = [1e-6]
    _max_value_list: List[float] = [1e6]
    
    
    def __init__(self, C_F: float = 300):
        """
        Equivalent circuit element: Ideal capacitor
        
        Model function for ideal capacitors
        Appearance: Vertical line (negative, up)
        Model:      Z = 1 / (jwC)

        Args:
            C_F (float): Capacitance in F
        """
        self.set_parameters(C_F)
        
    
    def evaluate(self, frequency_Hz: npt.ArrayLike) -> npt.ArrayLike:
        return 1 / (self.C_F * 1j * 2 * np.pi * frequency_Hz)
    
    
    @property 
    def C_F(self):
        return self._parameter_list[0]
    
    @C_F.setter
    def C_F(self, value):
        self._parameter_list[0] = value