import json
import numpy as np
import qubovert as qv
import matplotlib.pyplot as plt
import time as ooooo

with open('C:/Users/Михаил/Desktop/ChemKin/second_order.json') as f:
    data = json.load(f)

time = np.array(data['time'])
con = np.array(data['concentrations'])

N_add_k_p = []

N = len(time) - 1
M = len(con[0]) # 2
K = 1  # number of k's
n = 4  # digits in k  + 1
# dt = 0.5

N_add = 4 # number of additional dots between each pair + 1
        # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!
dT = 1

dt = dT/N_add

k_min = 0.125
k_max = 2 # here we should change


def recur(dt, c_AJ, c_BJ, N, N_add, ans, c_AJ_array, c_BJ_array):
    if N == 0:
        ans.append([c_AJ_array, c_BJ_array])
    else:
        recur(dt, c_AJ, c_BJ, N - 1, N_add, ans, c_AJ_array, c_BJ_array)
        c_A_new = ans[-1][0]
        c_B_new = ans[-1][1]
        for i in range(2**(N_add)):
            t1 = np.pad(ans[-1][0], (1+i, 0), 'constant', constant_values=(0,))[:-(i+1)] #тут ошибка.
            t2 = np.pad(ans[-1][0], (1+i, 0), 'constant', constant_values=(0,))[:-(i+1)]
            c_A_new = c_A_new - dt*ans[-1][1][i]*t1
            c_B_new = c_B_new - dt*ans[-1][1][i]*t2 #+here
        ans.append([c_A_new, c_B_new])



def c_wave(dt, c_AJ, c_BJ, N_add, k_vars, k_min, k_max, n):
    c_AJ_array = np.array([c_AJ] + [0]*(2**(N_add) - 1))
    c_BJ_array = np.array([c_BJ] + [0]*(2**(N_add) - 1))
    ans = []
    recur(dt, c_AJ, c_BJ, N_add+1, N_add, ans, c_AJ_array, c_BJ_array)
    
    # ((-1)**(j+1))

    c_A = ans[-1][0][0]
    c_B = ans[-1][1][0]
    for j in range(1,2**(N_add)):
        c_A += ans[-1][0][j]*(k_min**(j) + (k_max**(j) - k_min**(j))*sum((0.5 ** i) * k_vars[j-1][i] for i in range(1,n))) # тут знак
        c_B += ans[-1][0][j]*(k_min**(j) + (k_max**(j) - k_min**(j))*sum((0.5 ** i) * k_vars[j-1][i] for i in range(1,n)))
    return [c_A, c_B]


asdfasdf = ooooo.time()

k_vars = [[qv.PUBO.create_var('k'+str(i) + '_' + str(j)) for i in range(n)] for j in range(2**(N_add) - 1)] 
k_powers = [- (k_min**(j)+(k_max**(j)-k_min**(j))*sum((0.5 ** i) * k_vars[j-1][i] for i in range(1,n))) for j in range(1, 2**(N_add) - 1)]

# print('--------------------', k_powers)

# print('----------0000-------',(k_powers[2] - k_powers[0]**(2+1))**2)

F = 0
for i in range(1, N + 1):
    dt = (time[i]-time[i-1]) / N_add
    c_A = c_wave(dt, con[i - 1][0], con[i - 1][1], N_add, k_vars, k_min, k_max, n)[0]
    c_B = c_wave(dt, con[i - 1][0], con[i - 1][1], N_add, k_vars, k_min, k_max, n)[1]
    print(i, N) 
    F += (con[i][0] - c_A) ** 2 + (con[i][1] - c_B) ** 2

# reg_koeff = 123123

# F += reg_koeff*sum((k_powers[j] - k_powers[0]**(j+1))**2 for j in range(1, 2**(N_add) - 2))



int_F = F.to_qubo()
res = qv.sim.anneal_qubo(F, num_anneals=5000)

F_solution = res.best.state

k_p = [-i.value(F_solution) for i in k_powers]
print('N_add = ', N_add)
print('k_pows =',k_p)

# N_add_k_p.append([N_add, k_p])
# print('n = ', n)
# print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=', N_add_k_p)

# ppppppppppp = [val**(1/(i+1)) for i, val in enumerate(k_p)]

# print('data_k = ', ppppppppppp)

# print('k_aver',sum(ppppppppppp) / len(ppppppppppp))

# print('t =', ooooo.time() - asdfasdf)

# print(N_add_k_p)