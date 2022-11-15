# Import Packages
import numpy as np
import matplotlib.pyplot as plt
import math

# User Defined Functions

# Initialize
num_div = 30  # Number of divisions in x and y
# phi = 5.0*exp(-1500 * ((x - 0.25) ** 2 + y ** 2))  # Plotted equation
x = np.linspace(-1, 1, num_div)
y = np.linspace(-1, 1, num_div)
sol = np.zeros((num_div, num_div))

# Solve
for i in range(num_div - 1):
    for j in range(num_div - 1):
        sol[i][j] = 5.0 * math.exp(-1500 * ((x[i] - 0.25) ** 2 + y[j] ** 2))
