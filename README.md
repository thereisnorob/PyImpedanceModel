# PyImpedanceModel
PyImpedanceModel is a lightweight python framework for modeling electrical impedance models.

![Example data fit](https://picr.eu/images/2023/03/11/C2kEe.png)

In its current form, it can be used to setup generic equivalent-circuit models (e.g. for generating dummy electrical impedance spectroscopy (EIS) data) 
and to fit models to real-world impedance data.

Simple example scripts for a) seeting up arbitrary models, and b) fit a model to real-world data can be found in the repo.

## Usage
In its most basic form, a model can be initialized as follows:
```python
from PyImpedanceModel import *

model = ImpedanceModel([R(), (R(), C())])
```

Model notation:
Model elements in square brackets will be interpreted as connected in series. Elements in round brackets will be treated as connected in parallel. Both types can be combined to set up more complex models.

Example for simple model containing a resistor in series with a capacitor (R-C):
```python
model = ImpedanceModel([R(), C()])
```

Example for a simple parallel resistor-capacitor circuit:
```python
model = ImpedanceModel((R(), C()))
```

A standard Randles model (https://en.wikipedia.org/wiki/Randles_circuit) can be set up as follows:
```python
model = ImpedanceModel([R(20), (C(0.000025), [R(100), W(300)])])
model.evaluate(frequency_Hz = ImpedanceModel.get_log_frequency_range(800, 1))
```
The evaluate method allows to pass a user-defined frequency. In this case, this is used to rebuild the Randles circuit from Wikipedia.

## Implemented model elements
The following equivalent circuit elements are added:

### Ideal resistor `R(R_Ohm = 0.01)`
An ideal Ohmic resistor with resistance R in Ohm. 

### Ideal capacitor `C(C_F = 300)`
An ideal capacitor with capacitance C in F. 

### Ideal inductor `L(L_H = 1e-12)`
An ideal inductor with inductance L in H. 

### Constant-phase element (CPE) for non-ideal capcitors `CPE(Q0_Ohm_p_s_n = 30, n_1 = 0.95)`
In most cases ideal capacitors are insufficient for describing real-world capacitive impedance behaviors.

For such cases, a constant-phase element (short CPE) can be implemented. A very good introduction to CPEs can be found here: [https://lithiuminventory.com/experimental-electrochemistry/eis/constant-phase-element/index.html]

### CPE in parallel to a resistor `R_CPE(R_Ohm: = 0.01, tau_s: = 1, n_1: = 0.95)`
As CPEs are often used in parallel with an Ohmic resistor (often called ZARC element), this element can be used to model such a circuit.
It also replaces the capacitance parameter with a time constant (tau = R*C) which really helps understanding the connection to distribution of relaxation times (DRT). 

A brief ressource on ZARC elements can be found here:
[https://www.zimmerpeacocktech.com/2022/02/18/the-zarc-element/]

### Semi-inifite diffusion - Warburg `W(Aw_Ohm_p_s_sqrt = 1e-6)`
Warburg elements are used to model semi-infinite, linear diffusion behaviors in impedance spectra.
Its basically a frequentcy-independent 45° CPE. Aw is the Warburg constant.

### Finite-length, closed, or transmissive Warburg `FLW(Aw_Ohm_p_s_sqrt: = 1e-7, B_1_p_s_sqrt: = 1e-7)`
This element can be used to model transmissive mass-transport through thin layers with a finite length
(see https://lithiuminventory.com/experimental-electrochemistry/eis/diffusion-impedance/index.html).

### Finite-space, open, or reflective Warburg `FLW(Aw_Ohm_p_s_sqrt: = 0.001, B_1_p_s_sqrt: = 20)`
This element can be used to model reflective diffusion as it appears in porous electrodes
(see https://lithiuminventory.com/experimental-electrochemistry/eis/diffusion-impedance/index.html).
