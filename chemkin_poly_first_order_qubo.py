import json
import numpy as np
import scipy.special
import qubovert as qv
import matplotlib.pyplot as plt
import time
from utils import get_mininization_func

start_time = time.time()

with open('first_order.json') as f:
    data = json.load(f)

exp_time = np.array(data['time'])
con = np.array(data['concentrations'])
N_exp = len(exp_time)
n = 9  # digits in k  + 1

N_add = 10 # number of additional dots between each pair of experimental dots
          # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!

reg_koeff = 100000

k_vars = [[qv.QUBO.create_var('k%d%d' %(i, j)) for i in range(n)] for j in range(N_add + 2)]
k_vars = np.array(k_vars)
print("k_vars_shape = ", k_vars.shape)
k_pows = list()
for j in range(N_add + 2):
    k_pows.append(((-1)**j)*sum((0.5 ** i) * k_vars[j][i] for i in range(n)))

F = get_mininization_func(con, exp_time, k_pows, N_exp, N_add)

F += reg_koeff*(k_pows[0] - 1)**2
F += reg_koeff*sum((k_pows[j] - k_pows[1])**2 for j in range(1, N_add+2))

int_F = F.to_qubo()
res = qv.sim.anneal_qubo(F, num_anneals=1000)
end_time = time.time()
print("--- %s seconds ---" % (end_time - start_time))
print(res.best)

F_solution = res.best.state
print("\nF_solution =", F_solution)
print("k = ", -k_pows[0].value(F_solution))
print("F(k) =", F.value(F_solution))

plt.hist([sample.value for sample in res])
plt.xlabel("$F$")
plt.ylabel("counts")
plt.title("Values of $F$ sampled")
plt.show()
