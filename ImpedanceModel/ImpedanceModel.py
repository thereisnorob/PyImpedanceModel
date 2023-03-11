from __future__ import annotations
from typing import List, Tuple, Iterable
from .ImpedanceModelElement import ImpedanceModelElement
from scipy import optimize
from scipy.optimize import OptimizeResult
from time import perf_counter
import numpy as np
import numpy.typing as npt
import warnings

# Ignore RuntimeWarnings for overflowing numbers while settig up frequencies
warnings.filterwarnings('ignore', category = RuntimeWarning)

class ImpedanceModel:
    """
    Class for setting-up impedance models and fit them to measured spectra.
    
    @version:   AB-20230311
                (AA-20221224)
    @author:    Robert Leonhardt <mail@robertleonhardt.de>
    """
    
    # List with all model elements
    # Structe contains lists and tuples of ImpedanceModelElement
    _model_structure:  ImpedanceModelElement | Iterable[ImpedanceModelElement]
    _flat_model_list:   [ImpedanceModelElement]
    
    
    def __init__(self, model: ImpedanceModelElement | Iterable[ImpedanceModelElement], **kwargs):
        """
        Method to setup an impedance model based on ImpedanceModelElement subclass objects.
        
        Args:
            model_structure (ImpedanceModelElement | Iterable[ImpedanceModelElement]): Model structure.
        """
        self.__dict__.update(kwargs)
        
        # Store model structure
        self._model_structure   = model
        
        # Setup default (fallback) frequency range for evaluation
        self._frequency_Hz      = ImpedanceModel.get_log_frequency_range(1_000_000, 0.001)
        
        # Convert model sturcture to a list if its only one element
        if isinstance(self._model_structure, ImpedanceModelElement):
            self._model_structure = [self._model_structure]
        
        # Flatten model so that the parameters can be accessed during fitting process
        self._flat_model_list   = list(self._recursive_element_flattening(self._model_structure))
        
        # Evaluate model right here (so that all parameters can be used)
        _                       = self.evaluate()
        
        
    def _recursive_element_flattening(self, model_structure: ImpedanceModelElement | Iterable[ImpedanceModelElement]) -> List[ImpedanceModelElement]:
        """
        Method to recursivley flatten the model structure. This enables object-access to fitting algorithm.

        Args:
            model_structure (ImpedanceModelElement | Iterable[ImpedanceModelElement]): Model structure as defined in initializer.

        Returns:
            List[ImpedanceModelElement]: Flattened list of element in structure
        """             
        # Iterate through elements and (if element is iterable) unpack it recursively
        # NOTE: https://stackoverflow.com/a/14491059/12136370
        for element in model_structure:
            try:
                yield from self._recursive_element_flattening(element)
            except TypeError:
                yield element
                
            
    def evaluate(self, frequency_Hz: npt.ArrayLike = None, 
                 store_resulting_impedance_data = True) -> (npt.ArrayLike, npt.ArrayLike, npt.ArrayLike):
        """
        Method to evaluate impedance model (-> get impedance data)

        Args:
            frequency_Hz (npt.ArrayLike): Frequency list for which the model shall be evaluated (default to None -> default will be used)
            store_resulting_impedance_data (bool): If True, calculated impedance will be stored in model_obj.z_Ohm (etc.)

        Returns:
            npt.ArrayLike: Frequency (the same as in the input, just as a pass through)
            npt.ArrayLike: Impedance 
            npt.ArrayLike: Phi -> Phase shift
        """
        
        # Define a frequency list if none is given
        if frequency_Hz is None:
            frequency_Hz = self._frequency_Hz
        
        # Setup lists with data and recursively iterate through model structure
        z_Ohm = np.zeros(shape = len(frequency_Hz), dtype = np.complex128)
        z_Ohm = self._recursive_element_evaluation(z_Ohm, frequency_Hz, self._model_structure)
        
        # Store
        if store_resulting_impedance_data:
            self._frequency_Hz  = frequency_Hz
            self._z_Ohm         = z_Ohm 
            self._phi_deg       = np.angle(z_Ohm, deg = True)
        
        return frequency_Hz, z_Ohm, np.angle(z_Ohm, deg = True)
    
    
    def evaluate_with_parameters(self, frequency_Hz: npt.ArrayLike, parameter_list: List[float],
                                 store_resulting_impedance_data = True) -> (npt.ArrayLike, npt.ArrayLike, npt.ArrayLike, npt.ArrayLike):
        """
        Method to evaluate the model with specified parameters.

        Args:
            frequency_Hz npt.ArrayLike: List with frequencies for which the model shall be evaluated
            parameter_list List[float]: List with parameters (flat list)

        Returns:
            npt.ArrayLike: Frequency (the same as in the input, just as a pass through)
            npt.ArrayLike: Z -> Calculated impedance 
            npt.ArrayLike: Phi -> Phase shift
        """        
        # Set parameters and store original as they will be rerolled later
        original_parameters = self._get_model_parameters()
        self._set_model_parameters(parameter_list)
        
        # Evaluate model
        evaluation_result   = self.evaluate(frequency_Hz, store_resulting_impedance_data = store_resulting_impedance_data)
        
        # Reset parameters 
        self._set_model_parameters(original_parameters)
        
        return evaluation_result
    
                
    def _recursive_element_evaluation(self, result: npt.ArrayLike, frequency_Hz: npt.ArrayLike, elements: ImpedanceModelElement | Iterable[ImpedanceModelElement]) -> npt.ArrayLike:
        """
        Method to recursively iterate through the given model sturcture to calculate resulting impedance
        NOTE: Internal function ...

        Args:
            result (npt.ArrayLike): Result that the currently calculated part-result shall be added to
            frequency_Hz (npt.ArrayLike): Frequency list that the elements need to be evaluated
            elements (ImpedanceModelElement | Iterable[ImpedanceModelElement]): Element structure

        Returns:
            npt.ArrayLike: Calculated impedance
        """
        
        # Single Element -> Evaluate single element
        if isinstance(elements, ImpedanceModelElement):
            return elements.evaluate(frequency_Hz)
        
        # List -> Elements are in series
        # For in-series impedances: Z_res = sum(Z)
        if isinstance(elements, list):
            return sum([self._recursive_element_evaluation(result, frequency_Hz, element) for element in elements])
        
        # Tuple -> Elements are in-parallel
        # For in-parallel impedances: Z_res = 1 / sum(1 / Z)
        if isinstance(elements, tuple):
            return 1 / sum([1 / self._recursive_element_evaluation(result, frequency_Hz, element) for element in elements])
        
        
    def _get_model_parameters(self) -> List[float]:
        """
        Method to get all model parameters in a flat list

        Returns:
            List[float]: List with model parameters
        """
        # Return flattened parameters lists for each equivalent circuit element         
        return [parameter for parameter_list in self._flat_model_list for parameter in parameter_list.get_parameters()]
    
    
    def _set_model_parameters(self, parameter_list: List[float]) -> None:
        """
        Method to set model element parameters based on a flat list.
        This is useful for fitting algorithms etc.

        Args:
            parameter_list (List[float]): List with parameters
        """        
        number_of_parameters = len(self._get_model_parameters())
        
        # Error, when number of given parameters is insufficient
        if len(parameter_list) != number_of_parameters:
            raise Exception(f'Insufficient number of given model parameters ({len(parameter_list)} given but {number_of_parameters} required)')
        
        index = 0
        
        # Iterate through elements
        for element in self._flat_model_list:
            number_of_parameters = len(element.get_parameters())
            
            # Set parameters
            element.set_parameters(*parameter_list[index:(index + number_of_parameters)])
            
            index += number_of_parameters
            
            
    def fit_data(self, frequency_data_Hz: npt.ArrayLike, z_data_Ohm: npt.ArrayLike, 
                 p_0: List[float] = None, p_min: List[float] = None, p_max: List[float] = None,
                 apply_result: bool = True) -> (OptimizeResult, float):
        """
        Method for fitting the model to measurement data

        Args:
            frequency_data_Hz (npt.ArrayLike): Frequency from measurement data in Hz
            z_data_Ohm (npt.ArrayLike): Impedance data from measurement data in Ohm
            p_0 (List[float], optional): Initial parameter guesses. Defaults to None.
            p_min (List[float], optional): Minimum parameter boundaries. Defaults to None.
            p_max (List[float], optional): Maximum parameter boundaries. Defaults to None.
            apply_result (bool, optional): If true, parameters will be applied and the model will be evaluated. Defaults to True.

        Returns:
            OptimizeResult: SciPy optimize result object with parameters and residuals
            float: Fit time in s
        """
                
        # Set initial guess if none are given
        if not p_0:
            p_0 = self._get_model_parameters()
            
        # Set boundaries if none are given
        if not p_min:
            p_min = [boundary for boundary_list in self._flat_model_list for boundary in boundary_list._min_value_list]
        if not p_max:
            p_max = [boundary for boundary_list in self._flat_model_list for boundary in boundary_list._max_value_list]
        
        t_before = perf_counter()
        
        # Fit model
        fit = optimize.least_squares(
            fun         = self._get_residuals, 
            x0          = p_0, 
            args        = (frequency_data_Hz, z_data_Ohm), 
            bounds      = (p_min, p_max), 
            max_nfev    = 1000, 
            xtol        = 1e-6, 
            ftol        = None, #1e-12,
            gtol        = None)
        
        t_after     = perf_counter()
        fit_time_s  = t_after - t_before
        
        # Apply fitted parameters
        if apply_result:
            self._set_model_parameters(fit.x)
            self.evaluate(frequency_data_Hz)
        
        return fit, fit_time_s
    
    
    def _get_residuals(self, parameter_list: List[float], frequency_data_Hz: npt.ArrayLike, z_data_Ohm: npt.ArrayLike) -> List[float]:
        """
        Residual function that is used by the curve-fitting algorithm

        Args:
            parameter_list (List[float]): Input parameters (estimated by curve-fitting algorithm)
            frequency_data_Hz (npt.ArrayLike): Frequency list of input data
            z_data_Ohm (npt.ArrayLike): Impedance input data

        Returns:
            List[float]: Residual list
        """
        
        # Run model
        _, z_model_Ohm, _   = self.evaluate_with_parameters(frequency_data_Hz, parameter_list, store_resulting_impedance_data = False)
        
        # Calculate difference
        diff                = np.array(z_model_Ohm - z_data_Ohm)
        
        # Split real and imaginary part and merge dem to 1d-vector
        z1d                 = np.zeros(z_data_Ohm.size * 2, dtype = np.float64)
        z1d[0:z1d.size:2]   = diff.real
        z1d[1:z1d.size:2]   = diff.imag
        
        # Return difference as list
        return z1d
    
    
    @staticmethod
    def get_log_frequency_range(f_max_Hz: float, f_min_Hz: float, points_per_decade: int = 10) -> List[float]:
        """
        Method to generate frequency ranges

        Args:
            f_max_Hz (float): Maximum frequency (fast processes, induction)
            f_min_Hz (float): Minimum frequency (slow processes, diffusion)
            points_per_decade (int, optional): Number of points per decade. Defaults to 10.

        Returns:
            List[float]: Frequency range from max to min.
        """
        
        # Calculate number of points
        number_of_points = int(-np.floor(np.log10(f_min_Hz)) + np.floor(np.log10(f_max_Hz))) * points_per_decade
        
        # Make log-spaced ranged an return
        return np.flipud(np.geomspace(f_min_Hz, f_max_Hz, num = number_of_points, dtype = np.float64))
    
    
    @property
    def frequency_Hz(self) -> npt.ArrayLike:
        return np.array(self._frequency_Hz)
    
    @property
    def z_Ohm(self) -> npt.ArrayLike:
        return np.array(self._z_Ohm)
    
    @property 
    def phi_deg(self) -> npt.ArrayLike:
        return np.array(self._phi_deg)