# Import Packages
import numpy as np
import matplotlib.pyplot as plt
# User Defined Functions

# Mech 587 Problem Set 2
# Christian Rowsell (40131393)

# # Start values
# tfinal = 1.0  # Final time
# tinitial = 0.0  # Initial time
# nstep = 100  # Number of timesteps
# len = 1.0  # Length in x
# numx = 20  # Number of divisions in x
#
#
# def burgers(numx, len, tfinal, tinitial, nstep):
#     # Initialize
#     dt = (tfinal - tinitial) / nstep  # Timestep
#     dx = len / numx  # Spatial step
#     u = np.zeros(numx)  # Unknown solution vector
#     u_temp = np.zeros(numx)  # Temporary solution vector
#     x = np.zeros(numx)  # Position Vector
#     t = 0  # Current time
#     nu = 0.01
#     c = 1.0
#     tsteps = np.zeros(numx)  # Time vector
#     # Initial conditions
#     for i in range(0, numx):
#         u[i] = 0.5*np.sin(np.pi*2*dx*i)  # Apply initial conditions
#         x[i] = i * dx  # Fill position vector
#     u_temp = u  # Store temporary solution
#     # Loop over time
#     for i in range(1, numx - 1):
#         u[i] = u_temp[i] - 0.5 * (dt / dx) * c * (u_temp[i + 1] - u_temp[i - 1]) + nu * (dt / dx ** 2) * (u_temp[i + 1] - 2 * u_temp[i] - u_temp[i - 1])
#     u[numx] = u_temp[numx] - 0.5 * (dt / dx) * c * (u_temp[1] - u_temp[numx-1]) + nu * (dt / dx ** 2) * (u_temp[1] - 2 * u_temp[numx] + u_temp[numx - 1])
#     u[0] = u[numx]
#
#
#     plt.plot(x, u)
#     plt.show()
#
#
# burgers(numx, len, tfinal, tinitial, nstep)
# Import packages
import numpy as np
import matplotlib.pyplot as plt
import os
# User defined packages

# Path to save files
current_path = os.getcwd()
plot_folder = current_path + '/plots'

def burgers(numx, dom_len, tfinal, tinitial, dt):
    dx = dom_len / (numx - 1)  # Spatial step size
    x = np.zeros(numx)  # Position vector
    t = np.arange(tinitial, tfinal + dt, dt)  # tfinal + dt is needed to include tfinal in interval
    u = np.zeros((numx, len(t)))  # Solution vector
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
        u[-1, n+1] = u[j, n + 1]

        for j in range(1, numx - 1):
            #  Solve solution vector
            u[j, n + 1] = (dt * nu / dx ** 2) * (u[j + 1, n] - 2 * u[j, n] + u[j - 1, n])\
                          - (dt / 2 * dx) * u[j, n] * (u[j + 1, n] - u[j - 1, n]) + u[j, n]
    plt.plot(x, u[:, -1], label='Diffusion-Advection', color='b')
    plt.grid()
    plt.legend(loc='best')
    plt.title("Question 3 Problem Set 2, dt=%g, Number of grid points=%i" % (dt, numx))
    plt.xlabel("Position X")
    plt.ylabel("Solution U")
    save_folder = plot_folder + '_Num_X=%i' % numx
    if not os.path.exists(save_folder):
        os.makedirs(save_folder, exist_ok=True)
    plt.savefig(save_folder + '/dt=%g_Numx=%i' % (dt, numx) + '.png')
    plt.close()
    # plt.show()

numx = [20, 40]  # Number of discretization points
dom_len = 1.0  # Domain size
# dt = np.linspace(5e-3, 5e-1, 20)
dt = np.arange(5e-3, 5e-1 + 5e-3, 5e-3)
# dt = 5e-3  # Time step
tfinal = 1.0
tinitial = 0.0
# burgers(20, 1, 1, 0, 0.5)
for t in dt:
    for num in numx:
        print(t, numx)
        burgers(num, dom_len, tfinal, tinitial, t)
