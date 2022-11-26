# Import Packages
import numpy as np
import matplotlib.pyplot as plt


# User Defined Functions

# Initialize
num_div = 100  # Number of divisions
x_list = np.linspace(-1.0, 1.0, num_div)  # x coordinates
y_list = np.linspace(-1.0, 1.0, num_div)  # y coordinates
# Polar coordinates
x, y = np.meshgrid(x_list, y_list)  # Mesh
# Calculate phi
phi = 5.0 * np.exp(-1500.0 * ((x - 0.25) ** 2 + y ** 2))  # Initial phi

# Magnitude
u = y
v = -x
r = np.sqrt(u ** 2 + v ** 2)
theta = np.arctan(v/u)
# plt.contourf(x_list, y_list, u)
# plt.contourf(x_list, y_list, v)
# plt.contourf(x_list, y_list, r)
plt.contourf(x_list, y_list, theta)
# Adding
velocity = -x + y
# plt.contourf(x_list, y_list, velocity)

# Show plot
plt.xlabel('X')
plt.ylabel('Y')
plt.title('theta Velocity')
plt.colorbar()
plt.show()
