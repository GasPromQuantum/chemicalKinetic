import json
import numpy as np
import qubovert as qv
import time

start_time = time.time()

with open('Data/ChemKin/second_order.json') as f:
    data = json.load(f)

exp_time = np.array(data['time'])
con = np.array(data['concentrations'])
N = len(exp_time) - 1
M = len(con[0]) # 2
K = 1  # number of k's
n = 5  # digits in k  + 1
dt = 0.5

N_add = 2 # number of additional dots between each pair + 1
          # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!
dT = 1
k_min = 1
k_max = 2


def recur(dt, c_AJ, c_BJ, N, N_add, ans, c_AJ_array, c_BJ_array):
    if N == 0:
        ans.append([c_AJ_array, c_BJ_array])
    else:
        recur(dt, c_AJ, c_BJ, N - 1, N_add, ans, c_AJ_array, c_BJ_array)
        c_A_new = c_AJ_array
        c_B_new = c_BJ_array
        for i in range(N_add + 1):
            t1 = np.pad(c_AJ_array, (1+i, 0), 'constant', constant_values=0)[:-(i+1)]
            t2 = np.pad(c_AJ_array, (1+i, 0), 'constant', constant_values=0)[:-(i+1)]
            c_A_new = c_A_new - dt*c_BJ_array[i]*t1
            c_B_new = c_B_new - dt*c_BJ_array[i]*t2
            # np.roll(c_AJ_array, 1+i)
        ans.append([c_A_new, c_B_new])


def c_wave(dt, c_AJ, c_BJ, N, N_add, k_vars, k_min, k_max, n):
    c_AJ_array = np.array([c_AJ] + [0]*N_add)
    c_BJ_array = np.array([c_BJ] + [0]*N_add)
    ans = []
    recur(dt, c_AJ, c_BJ, N, N_add, ans, c_AJ_array, c_BJ_array)

    c_A = ans[-1][0][0]
    c_B = ans[-1][1][0]
    for j in range(1, N_add+1):
        c_A += ans[-1][0][j]*(k_min + (k_max - k_min)*sum((0.5 ** i) * k_vars[j-1][i] for i in range(1, n))) # тут минус???? Откуда он был раньше?
        c_B += ans[-1][1][j]*(k_min + (k_max - k_min)*sum((0.5 ** i) * k_vars[j-1][i] for i in range(1, n)))
    return [c_A, c_B]

k_vars = [[qv.PUBO.create_var('k%d' % i) for i in range(n)] for j in range(N_add)]

print(c_wave(dt, 1, 2, N, N_add, k_vars, k_min, k_max, n))

F = 0
for i in range(1, N + 1):
    dt = (exp_time[i]-exp_time[i-1]) / N_add
    c_A = c_wave(dt, con[i - 1][0], con[i - 1][1], int(dT / dt), N_add, k_vars, k_min, k_max, n)[0]
    c_B = c_wave(dt, con[i - 1][0], con[i - 1][1], int(dT / dt), N_add, k_vars, k_min, k_max, n)[1]
    F += (con[i][0] - c_A) ** 2 + (con[i][1] - c_B) ** 2

int_F = F.to_qubo()
res = qv.sim.anneal_qubo(F, num_anneals=1000)
end_time = time.time()
print("--- %s seconds ---" % (end_time - start_time))
print(res.best)