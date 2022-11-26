# Import functions
import numpy as np
import matplotlib.pyplot as plt
import math

# User defined functions

# Plot
# Exact Solution
Phi_exact_YX = np.genfromtxt("Phi_Exact_Y=X.csv", dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))
Phi_exact_Y05 = np.genfromtxt("Phi_Exact_Y=0.5.csv", dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))
Phi_exact_X05 = np.genfromtxt("Phi_Exact_X=0.5.csv",  dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))

# Grid 17
Phi_17_YX = np.genfromtxt("Phi17_Y=X.csv", dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))
Phi_17_Y05 = np.genfromtxt("Phi17_Y=0.5.csv", dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))
Phi_17_X05 = np.genfromtxt("Phi17_X=0.5.csv",  dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))

# Grid 65
Phi_65_YX = np.genfromtxt("Phi65_Y=X.csv", dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))
Phi_65_Y05 = np.genfromtxt("Phi65_Y=0.5.csv", dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))
Phi_65_X05 = np.genfromtxt("Phi65_X=0.5.csv",  dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))

# Grid 257
Phi_257_YX = np.genfromtxt("Phi257_Y=X.csv", dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))
Phi_257_Y05 = np.genfromtxt("Phi257_Y=0.5.csv", dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))
Phi_257_X05 = np.genfromtxt("Phi257_X=0.5.csv",  dtype=float, delimiter=',', skip_header=1, usecols=(4, 5))

# Plot Y=X Numerical Vs. Exact Solution
plt.plot(Phi_exact_YX[:, 0], Phi_exact_YX[:, 1], label='Exact Solution')
plt.plot(Phi_17_YX[:, 0], Phi_17_YX[:, 1], label='17 Grid Points')
plt.plot(Phi_65_YX[:, 0], Phi_65_YX[:, 1], label='65 Grid Points')
plt.plot(Phi_257_YX[:, 0], Phi_257_YX[:, 1], label='257 Grid Points')
plt.title("Exact Vs. Numerical Solution Along Line Y=X")
plt.xlabel("Arc Length")
plt.ylabel("Phi")
plt.legend(loc='best')
plt.grid()
plt.show()

# Plot Y=0.5 Numerical Vs. Exact Solution
plt.plot(Phi_exact_Y05[:, 0], Phi_exact_Y05[:, 1], label='Exact Solution')
plt.plot(Phi_17_Y05[:, 0], Phi_17_Y05[:, 1], label='17 Grid Points')
plt.plot(Phi_65_Y05[:, 0], Phi_65_Y05[:, 1], label='65 Grid Points')
plt.plot(Phi_257_Y05[:, 0], Phi_257_Y05[:, 1], label='257 Grid Points')
plt.title("Exact Vs. Numerical Solution Along Line Y=0.5")
plt.xlabel("Arc Length")
plt.ylabel("Phi")
plt.legend(loc='best')
plt.grid()
plt.show()

# Plot X=0.5 Numerical Vs. Exact Solution
plt.plot(Phi_exact_X05[:, 0], Phi_exact_X05[:, 1], label='Exact Solution')
plt.plot(Phi_17_X05[:, 0], Phi_17_X05[:, 1], label='17 Grid Points')
plt.plot(Phi_65_X05[:, 0], Phi_65_X05[:, 1], label='65 Grid Points')
plt.plot(Phi_257_X05[:, 0], Phi_257_X05[:, 1], label='257 Grid Points')
plt.title("Exact Vs. Numerical Solution Along Line X=0.5")
plt.xlabel("Arc Length")
plt.ylabel("Phi")
plt.legend(loc='best')
plt.grid()
plt.show()
# 6.93254896246994e-02,
# 1/65
L2 = [6.56526365343552e-02, 6.93254896246994e-02, 6.75070459558353e-02]
grid = [1/17, 1/65, 1/257]
slope = (math.log(L2[-1]) - math.log(L2[0])) / (math.log(grid[-1]) - math.log(grid[0]))
print(slope)
plt.loglog(grid, L2, label='L2 Norm Vs. Grid Size')
plt.grid(which='both')
plt.legend(loc='best')
plt.title("L2 Norm Vs. Grid Size")
plt.xlabel('Grid Size')
plt.ylabel("L2 Norm")
plt.show()
