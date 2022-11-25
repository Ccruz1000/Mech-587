# Import Packages
import numpy as np
import matplotlib.pyplot as plt


# User Defined Functions

# Initialize
num_div = 100  # Number of divisions
x_list = np.linspace(-1.0, 1.0, num_div)  # x coordinates
y_list = np.linspace(-1.0, 1.0, num_div)  # y coordinates
x, y = np.meshgrid(x_list, y_list)  # Mesh
# Calculate phi
phi = 5.0 * np.exp(-1500.0 * ((x - 0.25) ** 2 + y ** 2))  # Initial phi
# plt.contourf(x_list, y_list, phi)

# Not sure what b wants, either do velocity = u + v or do velocity = sqrt(u^2 + v^2)
# ask

# Magnitude
u = y
v = -x
vel = np.sqrt(u ** 2 + v ** 2)
plt.contourf(x_list, y_list, vel)
# Adding
velocity = -x + y
# plt.contourf(x_list, y_list, velocity)

# Show plot
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Phi')
plt.show()
