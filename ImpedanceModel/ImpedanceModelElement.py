from typing import List, Any
import numpy as np 
import numpy.typing as npt

class ImpedanceModelElement:
    """
    Parent class for equivalent circuit parameters
    
    @version:   AB-20230311
                (AA-20221224)
    @author:    Robert Leonhardt <mail@robertleonhardt.de>
    """
        
    # List with element parameters
    _parameter_list: List[float]
    
    
    def set_parameters(self, *parameters: List[float]) -> None:
        """
        Method to store parameters.
        
        Args:
            *parameters List[float]: List with parameters (ordered)
        """
        self._parameter_list = list(parameters)
        
    
    def get_parameters(self) -> List[float]:
        """
        Method to get parameters.

        Returns:
            List[float]: Parameter list
        """
        return self._parameter_list
    
    
    def evaluate(self) -> npt.ArrayLike:
        raise NotImplementedError()