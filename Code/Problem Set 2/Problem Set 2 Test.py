# Import packages
import numpy as np
import matplotlib.pyplot as plt
# User defined packages

numx = 40  # Number of discretization points
dom_len = 1.0  # Domain size
dx = dom_len / (numx - 1)  # Spatial step size
x = np.zeros(numx)  # Position vector
dt = 5e-3  # Time step
tfinal = 1.0
tinitial = 0.0
t = np.arange(tinitial, tfinal + dt, dt)  # tfinal + dt is needed to include tfinal in interval
u = np.zeros((numx, len(t)))  # Solution vector
u_new = np.zeros((numx, len(t)))  # Temporary solution vector
nu = 0.01
# Initial conditions
for i in range(numx):
    u[i, 0] = 0.5 * np.sin(np.pi * 2 * dx * i)  # Apply initial conditions
    x[i] = i * dx  # Fill position vector
for n in range(len(t) - 1):
    j = 0
    #  Solve initial point for periodic boundary condition
    u[j, n + 1] = (dt * nu / dx ** 2) * (u[j + 1, n] - 2 * u[j, n] + u[-1, n]) \
                  - (dt / 2 * dx) * u[j, n] * (u[j + 1, n] - u[-1, n]) + u[j, n]

    for j in range(1, numx - 1):
        #  Solve solution vector
        u[j, n + 1] = (dt * nu / dx ** 2) * (u[j + 1, n] - 2 * u[j, n] + u[j - 1, n])\
                      - (dt / 2 * dx) * u[j, n] * (u[j + 1, n] - u[j - 1, n]) + u[j, n]

plt.plot(x, u[:, -1])
plt.show()

