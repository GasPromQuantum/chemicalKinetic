import scipy
import scipy.integrate
import json
import numpy as np


y_0 = [0.93, 1.12, 2.56]

k1 = -1.1
k2 = -0.57

# для 3 ордер - в беседе

def fun(t, y):
    return [k1*y[0]*y[1], k1*y[0]*y[1]+k2*y[1]*y[2], k2*y[1]*y[2]]

sol = scipy.integrate.solve_ivp(fun, [0, 10], y_0, max_step = 0.1)

# print(sol.t, sol.y)

print(sol.t[::10])
# print(sol.y[0][::10], )

times = [i for i in sol.t[::10]]
cons = [[i, j, k] for i,j,k in zip(sol.y[0][::10], sol.y[1][::10], sol.y[2][::10])]

print(times)
print(cons)

data = {"time": times, "concentrations": cons}

json_object = json.dumps(data, indent=4)


with open('C:\\py\\QuantEuler\\third_order1.json', "w") as f:
    f.write(json_object)


import matplotlib.pyplot as plt

plt.plot(sol.t, sol.y[0], color='r', label='cA')
plt.plot(sol.t, sol.y[1], color='r', label='cB') 
plt.plot(sol.t, sol.y[2], color='r', label='cC')  
plt.legend()
plt.xlabel("T")
plt.ylabel("$C$")
# plt.title("Values of $F$ sampled")
plt.show()