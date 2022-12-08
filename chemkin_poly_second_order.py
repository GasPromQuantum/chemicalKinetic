import json
import numpy as np
import qubovert as qv
import matplotlib.pyplot as plt

with open('Data/ChemKin/second_order.json') as f:
    data = json.load(f)

time = np.array(data['time'])
con = np.array(data['concentrations'])
N = len(time) - 1
M = len(con[0])
K = 1  # number of k's
n = 5  # digits in k  + 1
dt = 0.5

N_add = 2 # number of additional dots between each pair + 1
          # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!
dT = 1
k_min = 0.125
k_max = 2

def recur(dt, c_AJ, c_BJ, N, k, ans):
    if N == 0:
        ans.append([c_AJ, c_BJ])
    else:
        recur(dt, c_AJ, c_BJ, N - 1, k, ans)
        temp = [ans[-1][0] - k * dt * ans[-1][0] * ans[-1][1], ans[-1][1] - k * dt * ans[-1][0] * ans[-1][1]]
        ans.append(temp)

def c_wave(dt, c_AJ, c_BJ, N, k):
    ans = []
    recur(dt, c_AJ, c_BJ, N, k, ans)
    return ans[-1]

k_vars = [qv.PUBO.create_var('k%d' % i) for i in range(n)]
k = - (k_min + (k_max - k_min)*sum((0.5 ** i) * k_vars[i] for i in range(1, n)))
print(k_vars)
print(k)

F = 0
for i in range(1, N + 1):
    dt = (time[i]-time[i-1]) / N_add
    c_A = c_wave(dt, con[i - 1][0], con[i - 1][1], int(dT / dt), k)[0]
    c_B = c_wave(dt, con[i - 1][0], con[i - 1][1], int(dT / dt), k)[1]
    F += (con[i][0] - c_A) ** 2 + (con[i][1] - c_B) ** 2

int_F = F.to_pubo()
res = qv.sim.anneal_pubo(F, num_anneals=1000)
print(res.best)