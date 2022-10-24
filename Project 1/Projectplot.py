import matplotlib.pyplot as plt
import math

L2 = [6.048073324687e-04, 1.561666503124e-04, 3.966651043739e-05]
grid = [1/16, 1/32, 1/64]

slope = (math.log(L2[-1]) - math.log(L2[0])) / (math.log(grid[-1]) - math.log(grid[0]))
print(slope)
plt.loglog(grid, L2, label='L2 Norm Vs. Grid Size')
plt.grid(which='both')
plt.legend(loc='best')
plt.title("L2 Norm Vs. Grid Size")
plt.xlabel('Grid Size')
plt.ylabel("L2 Norm")
plt.show()