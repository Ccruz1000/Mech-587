import numpy as np
import matplotlib.pyplot as plt

# Mech 587 - CFD
# Sample code for three integration methods

nstep = 2  # Number of time steps
dt = 1/nstep # Time step
# Initialize arrays
u1 = np.zeros(nstep)

u2 = np.zeros(nstep)
u3 = np.zeros(nstep)
uex = np.zeros(nstep)
t = np.zeros(nstep)
# Initial conditions
t[0] = 0
u1[0] = 1
u2[0] = 1
u3[0] = 1
uex[0] = 1
# Solve equations
for i in range(1, nstep):
    u1[i] = u1[i-1] - dt*u1[i-1]  # Forward Euler
    u2[i] = u2[i-1] / (1.0 + dt)  # Backward Euler
    u3[i] = u3[i-1] * (1.0 - (dt / 2)) / (1.0 + (dt/2))  # Trapezoidal Rule
    t[i] = t[i-1] + dt
    uex[i] = np.exp(-t[i])
#
plt.plot(t, u1, label='Forward Euler', color='blue')
plt.plot(t, u2, label='Backward Euler', color='red')
plt.plot(t, u3, label='Trapezoidal Rule', color='orange')
plt.plot(t, uex, label='Exact Solution', linestyle='--', color='black')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Solution')
plt.title('ODE Integrator Method Comparison')
plt.show()
