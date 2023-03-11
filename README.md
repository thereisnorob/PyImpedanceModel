# PyImpedanceModel
PyImpedanceModel is a lightweight python framework for modeling electrical impedance models.

In its current form, it can be used to setup generic equivalent-circuit models (e.g. for generating dummy electrical impedance spectroscopy (EIS) data) 
and to fit models to real-world impedance data.

Simple example scripts for a) seeting up arbitrary models, and b) fit a model to real-world data can be found in the repo.

## Usage:
In its most basic form, a model can be initialized as follows:
```
from ImpedanceModel import *

model = ImpedanceModel([R(), (R(), C())])
```

Model notation:
Model elements in square brackets will be interpreted as connected in series. Elements in round brackets will be treated as connected in parallel. Both types can be combined to set up more complex models.

Example for simple model containing a resistor in series with a capacitor (R-C):
```
model = ImpedanceModel([R(), C()])
```

Example for a simple parallel resistor-capacitor circuit:
```
model = ImpedanceModel((R(), C()))
```

A standard Randles model (https://en.wikipedia.org/wiki/Randles_circuit) can be set up as follows:
```
model = ImpedanceModel([R(20), (C(0.000025), [R(100), W(300)])])
model.evaluate(frequency_Hz = ImpedanceModel.get_log_frequency_range(800, 1))
```
The evaluate method allows to pass a user-defined frequency. In this case, this is used to rebuild the Randles circuit from Wikipedia.

## Implemented model elements
The following equivalent circuit elements are added:

### Ideal resistor `R()`
An ideal resistor with resistance R. 