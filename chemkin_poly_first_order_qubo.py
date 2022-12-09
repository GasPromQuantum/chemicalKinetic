import json
import numpy as np
import scipy.special
import qubovert as qv
import matplotlib.pyplot as plt
import time

start_time = time.time()

with open('Data/ChemKin/first_order.json') as f:
    data = json.load(f)

exp_time = np.array(data['time'])
con = np.array(data['concentrations'])
N = len(exp_time) - 1
M = len(con[0])
K = 1  # number of k's
n = 5  # digits in k  + 1

N_add = 2 # number of segments between [\tau_{i}, \tau_{i+1}]
          # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!
# dT = 1
# r = int(dT/dt)
k_min = 0.125
k_max = 2

k_vars = [[qv.PUBO.create_var("k" + str(i) + str(j)) for i in range(n)] for j in range(N_add + 1)]
k = [0] * (N_add + 1)
for j in range(N_add + 1):
    k[j] = - (k_min + (k_max - k_min)*sum((0.5 ** i) * k_vars[j][i] for i in range(1, n)))

F = 0
for i in range(1, N + 1):
    delta = con[i][0]
    dt = (exp_time[i]-exp_time[i-1]) / N_add
    for j in range(N_add + 1):
        C = scipy.special.binom(N_add + 1, j)
        delta -= C * con[i - 1][0] * (dt**j) * k[j]
    F += delta ** 2

int_F = F.to_qubo()
res = qv.sim.anneal_qubo(F, num_anneals=1000)
end_time = time.time()
print("--- %s seconds ---" % (end_time - start_time))
print(res.best)

F_solution = res.best.state
print("\nF_solution =", F_solution)
print("k =", k[0].value(F_solution))
print("F(k) =", F.value(F_solution))

plt.hist([sample.value for sample in res])
plt.xlabel("$F$")
plt.ylabel("counts")
plt.title("Values of $F$ sampled")
plt.show()
