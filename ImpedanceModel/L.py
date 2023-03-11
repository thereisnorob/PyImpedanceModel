from __future__ import annotations
import numpy as np
import numpy.typing as npt
from .ImpedanceModelElement import ImpedanceModelElement

class L(ImpedanceModelElement):
    """
    Equivalent circuit element: Ideal inductor
    
    @version:   AB-20230311
                (AA-20221224)
    @author:    Robert Leonhardt <mail@robertleonhardt.de>
    """
    
    # Default min/max values
    _min_value_list: List[float] = [1e-12]
    _max_value_list: List[float] = [1e4]
    
    
    def __init__(self, L_H: float = 1e-12):
        """
        Equivalent circuit element: Ideal inductor
        
        Model function for ideal inductors 
        Appearance: Vertical line (positive, down)
        Model:      Z = jwL

        Args:
            L_H (float): Inductance in H
        """
        self.set_parameters(L_H)
        
    
    def evaluate(self, frequency_Hz: npt.ArrayLike) -> npt.ArrayLike:
        return 1j * 2 * np.pi * frequency_Hz * self.L_H
    
    
    @property 
    def L_H(self):
        return self._parameter_list[0]
    
    @L_H.setter
    def L_H(self, value):
        self._parameter_list[0] = value