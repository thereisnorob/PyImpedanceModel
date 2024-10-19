import matplotlib.pyplot as plt
import numpy as np
from PyImpedanceModel import *

# Setup test models
model = ImpedanceModel(
    [R(0), 
     R_CPE(R_Ohm = 1, tau_s = 0.01, n_1 = 0.9), 
     R_CPE(R_Ohm = 8, tau_s = 1, n_1 = 1), 
     FSW(0.5, 20)]
)

# Setup plot
fig         = plt.figure(figsize = (9,4))
nyquist_ax  = fig.add_subplot()

# Plot data
nyquist_ax.plot(model.z_Ohm.real, -model.z_Ohm.imag, '-o', color = 'k', markerfacecolor = 'None', label = 'Model')

# Config plot
for ax in [nyquist_ax]:
    ax.axhline(0, color = 'k', lw = 0.5)
    ax.set_title('Impedance data')
    ax.set_xlabel('Re(z) in Ω')
    ax.set_ylabel('-Im(z) in Ω')
    ax.minorticks_on()
    ax.grid(visible = True, which = 'major', alpha = 0.4)
    ax.grid(visible = True, which = 'minor', alpha = 0.1)
    ax.set_aspect('equal')
    ax.legend()

plt.show()