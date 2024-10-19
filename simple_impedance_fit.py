import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PyImpedanceModel import *

# Setup test models
# As this model will be fitted to actual data, we won't bother coming up with initial values
model = ImpedanceModel(
    [(L(), R()), R(), R_CPE(), R_CPE(), W()]
)

# Read in data
impedance_data          = pd.read_csv('./example_eis_data.csv', skiprows = 1)
impedance_data['z_Ohm'] = impedance_data.z_real_Ohm + 1j * impedance_data.z_imag_Ohm

# Fit model to data (depending on the model, this can take some seconds ...)
fit, fit_time_s         = model.fit_data(impedance_data.frequency_Hz, impedance_data.z_Ohm)

# Setup plot
fig         = plt.figure(figsize = (9,4))
nyquist_ax  = fig.add_subplot()

# Plot data
nyquist_ax.plot(1000 * np.real(impedance_data.z_Ohm), -1000 * np.imag(impedance_data.z_Ohm), '-o', color = 'k', markerfacecolor = 'None', label = 'Data')
nyquist_ax.plot(1000 * model.z_Ohm.real, -1000 * model.z_Ohm.imag, '-+', color = '#ef233c', markerfacecolor = 'None', label = 'Fit')

# Config plot
for ax in [nyquist_ax]:
    ax.axhline(0, color = 'k', lw = 0.5)
    ax.set_title('Impedance data')
    ax.set_xlabel('Re(z) in mΩ')
    ax.set_ylabel('-Im(z) in mΩ')
    ax.minorticks_on()
    ax.grid(visible = True, which = 'major', alpha = 0.4)
    ax.grid(visible = True, which = 'minor', alpha = 0.1)
    ax.set_aspect('equal')
    ax.legend()

plt.show()