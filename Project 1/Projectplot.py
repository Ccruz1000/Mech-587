import matplotlib.pyplot as plt

L2 = [6.048073324687e-04, 1.561666503124e-04, 3.966651043739e-05]
grid = [289, 1089, 4225]

plt.loglog(grid, L2, label='L2 Norm Vs. Grid Size')
plt.grid()
plt.legend(loc='best')
plt.title("L2 Norm Vs. Grid Size")
plt.xlabel('Grid Size')
plt.ylabel("L2 Norm")
plt.show()