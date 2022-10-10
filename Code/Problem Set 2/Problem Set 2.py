# Import Packages
import numpy as np
import matplotlib.pyplot as plt
# User Defined Functions

# Mech 587 Problem Set 2
# Christian Rowsell (40131393)

# Start values
tfinal = 1.0  # Final time
tinitial = 0.0  # Initial time
nstep = 100  # Number of timesteps
len = 1.0  # Length in x
numx = 20  # Number of divisions in x


def burgers(numx, len, tfinal, tinitial, nstep):
    # Initialize
    dt = (tfinal - tinitial) / nstep  # Timestep
    dx = len / numx  # Spatial step
    u = np.zeros(numx)  # Unknown solution vector
    u_temp = np.zeros(numx)  # Temporary solution vector
    x = np.zeros(numx)  # Position Vector
    # Initial conditions
    for i in range(0, numx):
        u[i] = 0.5*np.sin(np.pi*2*dx*i)  # Apply initial conditions
        x[i] = i * dx  # Fill position vector
    plt.plot(x, u)
    plt.show()


burgers(numx, len, tfinal, tinitial, nstep)
