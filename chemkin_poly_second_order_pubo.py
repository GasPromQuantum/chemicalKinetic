import json
import numpy as np
import qubovert as qv
import time
import matplotlib.pyplot as plt

#start_time = time.time()

with open('second_order.json') as f:
    data = json.load(f)

exp_time = np.array(data['time'])
con = np.array(data['concentrations'])
N = len(exp_time) - 1
M = len(con[0])
#K = 1  # number of k's
n = 9  # digits in k  + 1
#dt = 0.5
anneals=1000

#dT = 1
k_min = 0.125
k_max = 2


def c_wave(dt, c_AJ, c_BJ, N, k):
    temp = [c_AJ, c_BJ]
    for _ in range(N):
        temp = [temp[0] + k * dt * temp[0] * temp[1], temp[1] + k * dt * temp[0] * temp[1]]

    return temp


k_vars = [qv.PUBO.create_var('k%d' % i) for i in range(n)]
k = - (k_min + (k_max - k_min)*sum((0.5 ** i) * k_vars[i] for i in range(1, n)))
print(k_vars)
print(k)
#N_add = 2 # number of additional dots between each pair + 1
          # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!
kF = []
#print(exp_time[2]-exp_time[1])
for N_add in range(7, 8):
    start_time = time.time()
    F = 0
    for i in range(1, N + 1):
        dt = (exp_time[i]-exp_time[i-1]) / N_add
        c_pair = c_wave(dt, con[i - 1][0], con[i - 1][1], N_add+1, k)
        c_A = c_pair[0]
        c_B = c_pair[1]
        F += (con[i][0] - c_A) ** 2 + (con[i][1] - c_B) ** 2
    int_F = F.to_pubo()
    res = qv.sim.anneal_pubo(F, num_anneals=anneals)
    end_time = time.time()
    print("--- %s seconds ---" % (end_time - start_time))
    print(res.best)
    F_solution = res.best.state
    print("k=",k.value(F_solution))
    kF.append([N_add,k.value(F_solution), F.value(F_solution)])

    plt.hist([sample.value for sample in res])
    plt.xlabel("$F$")
    plt.ylabel("counts")
    plt.title("Values of $F$ sampled")
    #plt.show()
    plt.savefig('second_pubo_n_%d_anneals_%d_N_add_%d.png' %(n,anneals,N_add))


    
print(kF)
