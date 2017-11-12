import matplotlib.pyplot as plt
import numpy as np
import photovoltaic as pv

IL = 0.035
I0 = 1e-12
V = np.linspace(0, 0.7, 100)
I = pv.I_cell(V, IL, I0)
plt.ylim(0, 0.04)
plt.xlabel('voltage (V)')
plt.ylabel('current (A)')
plt.grid(True)
plt.plot(V, I, label="Ideal Cell")
plt.legend(loc='lower left')
plt.savefig('plot_ideal_cell.png')
plt.show()
