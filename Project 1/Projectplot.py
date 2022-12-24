import matplotlib.pyplot as plt
import math

L2X = [2.558504873729e-02, 6.321108496806e-03, 1.569247085119e-03]
L2Y = [2.55850487372923e-02, 6.32110849680624e-03, 1.56924708511860e-03]
grid = [1/16, 1/32, 1/64]

slopex = (math.log(L2X[-1]) - math.log(L2X[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopey = (math.log(L2Y[-1]) - math.log(L2Y[0])) / (math.log(grid[-1]) - math.log(grid[0]))
print("X Error is equal to " + str(slopex))
print("Y Error is equal to " + str(slopey))
plt.loglog(grid, L2X, label='L2X Norm Vs. Grid Size')
plt.loglog(grid, L2Y, label='L2Y Norm Vs. Grid Size')
plt.grid(which='both')
plt.legend(loc='best')
plt.title("L2 Norm Vs. Grid Size")
plt.xlabel('Grid Size')
plt.ylabel("L2 Norm")
plt.show()
