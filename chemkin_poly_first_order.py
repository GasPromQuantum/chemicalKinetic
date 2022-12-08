import json
import numpy as np
import scipy.special
import qubovert as qv
import matplotlib.pyplot as plt

with open('Data/ChemKin/first_order.json') as f:
    data = json.load(f)

time = np.array(data['time'])
con = np.array(data['concentrations'])
N = len(time) - 1
M = len(con[0])
K = 1  # number of k's
n = 5  # digits in k  + 1

N_add = 2 # number of additional dots between each pair + 1
          # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!
# dT = 1
# r = int(dT/dt)
k_min = 0.125
k_max = 2

k_vars = [qv.PUBO.create_var('k%d' % i) for i in range(n)] #(1,n)???
k = - (k_min + (k_max - k_min)*sum((0.5 ** i) * k_vars[i] for i in range(1,n))) #(1,n)???
print(k_vars)
print(k)

temp = 1
k_pows = []
for i in range(N_add+1):
    k_pows.append(temp)
    temp *= k

F = 0
for i in range(1, N + 1):
    delta = con[i][0]
    dt = (time[i]-time[i-1]) / N_add
    for j in range(N_add + 1):
        C = scipy.special.binom(N_add, j)
        delta -= C * con[i - 1][0] * (dt**j) * k_pows[j]
    F += delta ** 2

print("F =", F)

int_F = F.to_pubo()
print("\nF vars:", F.variables)
print("int_F vars:", int_F.variables)
print("degree =",F.degree)

print("\nanneal_pubo solution:")
res = qv.sim.anneal_pubo(F, num_anneals=1000)
print(res.best)

F_solution = res.best.state
print("\nF_solution =", F_solution)
print("k =", k.value(F_solution))
print("F(k) =", F.value(F_solution))

plt.hist([sample.value for sample in res])
plt.xlabel("$F$")
plt.ylabel("counts")
plt.title("Values of $F$ sampled")
plt.show()