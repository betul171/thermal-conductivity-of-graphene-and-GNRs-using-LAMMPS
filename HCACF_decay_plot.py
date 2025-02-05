###### Make sure HCACF decays to 0.   ##############

import numpy as np
import matplotlib.pyplot as plt

# Load the flux correlation data from the file
data = np.loadtxt('a.dat', comments='#')  # Skip comments in the file if any
time = data[:, 1] # in ps
time = time * 0.001 # in fs
flux_corr_x = data[:, 3]  
flux_corr_y = data[:, 4]  
flux_corr_z = data[:, 5]  

plt.figure(figsize=(10, 6))

plt.plot(time, flux_corr_y, label='Flux Correlation in Y', linestyle='-')


plt.xlabel('Time (fs)')
plt.ylabel('Heat Current Auto-Correlation Function')
plt.title('HCACF vs Time for Different Directions')
plt.legend()
plt.grid(True)


plt.savefig('HCACF_vs_time.png') 
plt.show()  
