# Import Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
# Import User Defined Functions

# Mech 587 - CFD
# Christian Rowsell (40131393)

# Find path to save files
current_path = os.getcwd()
plot_folder = current_path + '/' + 'Part_A'

start_t = 0  # first t value
final_t = 8  # Last t value used


# Define function
def q3(dt, tinitial, tfinal):
    # Initialize arrays
    nstep = (tfinal - tinitial) / dt
    u1 = np.zeros(int(nstep))
    u2 = np.zeros(int(nstep))
    u3 = np.zeros(int(nstep))
    e1 = np.zeros(int(nstep))
    e2 = np.zeros(int(nstep))
    e3 = np.zeros(int(nstep))
    uex = np.zeros(int(nstep))
    t = np.zeros(int(nstep))

    # Initial conditions
    t[0] = tinitial
    u1[0] = 1
    u2[0] = 1
    u3[0] = 1
    uex[0] = 1
    # Solve equations
    for i in range(1, int(nstep)):
        u1[i] = u1[i-1] - 2*dt*u1[i-1]  # Forward Euler
        u2[i] = u2[i-1] / (1.0 + 2*dt)  # Backward Euler
        u3[i] = u3[i-1] * (1.0 - dt) / (1.0 + dt)  # Trapezoidal Rule
        t[i] = t[i-1] + dt
        uex[i] = np.exp(-2*t[i])

        # If code below works
        # if i == 4:
        #     err1 = abs(uex[i] - u1[i])  # Forward Euler
        #     err2 = abs(uex[i] - u2[i])  # Backward Euler
        #     err3 = abs(uex[i] - u3[i])  # Trapezoidal Rule
        #     return [err1, err2, err3]

    #     # Error Analysis
        e1[i] = abs(uex[i] - u1[i])  # Forward Euler
        e2[i] = abs(uex[i] - u2[i])  # Backward Euler
        e3[i] = abs(uex[i] - u3[i])  # Trapezoidal Rule
    #
    # # Interpolate equations to allow finding of error analysis at t=4 exactly
    eqn1 = interp1d(t, e1)
    eqn2 = interp1d(t, e2)
    eqn3 = interp1d(t, e3)
    # Finding Absolute Error
    err1 = eqn1(4)
    err2 = eqn2(4)
    err3 = eqn3(4)



    # Plot Equations
    # plt.plot(t, u1, label='Forward Euler', color='blue')
    # plt.plot(t, u2, label='Backward Euler', color='red')
    # plt.plot(t, u3, label='Trapezoidal Rule', color='orange')
    # plt.plot(t, uex, label='Exact Solution', linestyle='--', color='black')
    # plt.grid()
    # plt.legend(loc='best')
    # plt.xlabel('Time')
    # plt.ylabel('Solution')
    # title = 'Problem Set 1, Question 3, Part 1, dt = ' + str(dt)
    # plt.title(title)
    # # Save Files
    # if not os.path.exists(plot_folder):
    #    os.makedirs(plot_folder, exist_ok=True)
    # plt.savefig(plot_folder + '/' + title + '.png', bbox_extra_artists='legend_outside')
    # plt.show()
    # plt.close()
    # Dont return if code below works
    return err1, err2, err3

# Part 1
# dt_list = [0.1, 0.2, 0.4, 0.8]
# for dt in dt_list:
#     q3(dt, 0, 8)

# Part 2
# Code to find values for dt that will always pass 4

'''
DIDNT WORK
# value = []  # Placedholder
# dt_list = []  # Values that will pass 4
# increment = 0.0001
# end = finalt / increment
# for i in range(1, int(end)):
#     value.append(i * increment)  # Create list of all possible time steps based on increment
# for j in range(len(value)):
#     if (4 / value[j]) % 1 == 0:
#         dt_list.append(value[j])
#        If 4 divided by given time step is integer, 4 will be passed
#        and as such is stored in the array of suitable time steps
# err = []
for i in range(len(dt_list)):
    err.append(q3(dt_list[i], start_t, final_t))
'''

dt_list = np.linspace(0.001, 1, 10000)  # Take 10000 equally spaced divisions between 0 and 1
err1 = np.zeros(len(dt_list))
err2 = np.zeros(len(dt_list))
err3 = np.zeros(len(dt_list))

for i in range(len(dt_list)):
    err1[i], err2[i], err3[i] = q3(dt_list[i], start_t, final_t)

# plt.plot(dt_list, err1, label='Forward Euler', color='blue')
plt.plot(dt_list, err2, label='Backward Euler', color='red')
plt.plot(dt_list, err3, label='Trapezoidal Rule', color='orange')
plt.grid()
plt.title('Absolute Error Interpolated Value Comparison at t = 4.0')
plt.legend(loc='best')
plt.xlabel('Timestep')
plt.ylabel('Absolute Error')
plt.yscale('log')
plt.xscale('log')
plt.show()
