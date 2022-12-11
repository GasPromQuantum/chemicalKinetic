import json
import numpy as np
import qubovert as qv
import matplotlib.pyplot as plt
import time as ooooo


"""
eq-n:

c1' = k1 c1 c2
c2' = k1 c1 c2 + k2 c2 c3
c3' = k2 c2 c3
"""

with open('C:/py/QuantEuler/third_order1.json') as f:
    data = json.load(f)

time = np.array(data['time'])
con = np.array(data['concentrations'])


N = len(time) - 1
M = len(con[0]) # 2
K = 2  # number of k's
n = 3  # digits in k  + 1

N_add = 4 # number of additional dots between each pair + 1
        # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!
dT = 1

dt = dT/N_add


k_min = 0.125
k_max = 2


def shift(arr, i, j): #i раз вниз, j раз вправо
    temp = np.pad(arr, [(i, 0), (j, 0)], 'constant', constant_values=[(0,), (0,)])
    return temp[:len(temp) - i, :len(temp[0]) - j]



def recur(dt, c_AJ, c_BJ, c_CJ, N, N_add, ans, c_AJ_array, c_BJ_array, c_CJ_array):
    if N == 1: # почему-то он делает на 1 итерацию больше
        ans.append([c_AJ_array, c_BJ_array, c_CJ_array])
    else:
        recur(dt, c_AJ, c_BJ, c_CJ, N - 1, N_add, ans, c_AJ_array, c_BJ_array, c_CJ_array)
        c_A_new = ans[-1][0]
        c_B_new = ans[-1][1]
        c_C_new = ans[-1][2]

        for i in range(2**N_add):
            for j in range(2**N_add):
                c_A_new = c_A_new - ans[-1][1][i][j]*dt*(shift(ans[-1][0], i+1, j))# k1 c1(0) c2(1)        
                c_C_new = c_C_new - ans[-1][1][i][j]*dt*(shift(ans[-1][2], i, j+1))#k2 c2(1) c3(2)
                c_B_new = c_B_new - ans[-1][1][i][j]*dt*(shift(ans[-1][2], i, j+1))#k2 c2(1) c3(2)
                c_B_new = c_B_new - ans[-1][1][i][j]*dt*(shift(ans[-1][0], i+1, j))# k1 c1(0) c2(1)
        ans.append([c_A_new, c_B_new, c_C_new])

def c_wave(dt, c_AJ, c_BJ, c_CJ, N_add, k_vars, k_min, k_max, n): #если что дописать.
    c_AJ_array = np.zeros((2**(N_add), 2**(N_add)))
    c_BJ_array = np.zeros((2**(N_add), 2**(N_add)))
    c_CJ_array = np.zeros((2**(N_add), 2**(N_add)))

    c_AJ_array[0][0] = c_AJ
    c_BJ_array[0][0] = c_BJ
    c_CJ_array[0][0] = c_CJ

    ans = []
    recur(dt, c_AJ, c_BJ, c_CJ, N_add+1, N_add, ans, c_AJ_array, c_BJ_array, c_CJ_array)


    c_A = ans[-1][0][0][0]
    c_B = ans[-1][1][0][0]
    c_C = ans[-1][2][0][0]


    for i in range(1, 2**N_add):
        j = 0 

        c_A += ans[-1][0][j][i]*(k_min**(i+j)+(k_max**(i+j)-k_min**(i+j))*sum((0.5 ** ass) * k_vars[j][i][ass] for ass in range(1, n)))
        c_B += ans[-1][1][j][i]*(k_min**(i+j)+(k_max**(i+j)-k_min**(i+j))*sum((0.5 ** ass) * k_vars[j][i][ass] for ass in range(1, n)))
        c_C += ans[-1][2][j][i]*(k_min**(i+j)+(k_max**(i+j)-k_min**(i+j))*sum((0.5 ** ass) * k_vars[j][i][ass] for ass in range(1, n)))


    for j in range(1, 2**N_add):
        for i in range(2**N_add - j):
            c_A += ans[-1][0][j][i]*(k_min**(i+j)+(k_max**(i+j)-k_min**(i+j))*sum((0.5 ** ass) * k_vars[j][i][ass] for ass in range(1, n)))
            c_B += ans[-1][1][j][i]*(k_min**(i+j)+(k_max**(i+j)-k_min**(i+j))*sum((0.5 ** ass) * k_vars[j][i][ass] for ass in range(1, n)))
            c_C += ans[-1][2][j][i]*(k_min**(i+j)+(k_max**(i+j)-k_min**(i+j))*sum((0.5 ** ass) * k_vars[j][i][ass] for ass in range(1, n)))
    return [c_A, c_B, c_C]


k_vars = [[[qv.PUBO.create_var('k'+str(k)+'_'+str(j)+'_'+str(i)) for i in range(n)] for j in range((2**(N_add)) - k)] for k in range(2**(N_add))]

k = [[[(k_min**(j+k)+(k_max**(j+k)-k_min**(j+k))*sum((0.5 ** i) * k_vars[k][j][i] for i in range(1,n)))] for j in range((2**(N_add)) - k)] for k in range(2**(N_add))]


print(k)

# print('-----------', c_wave(dt, con[1 - 1][0], con[1 - 1][1], con[1-1][2], N_add, k_vars, k_min, k_max, n))


F = 0
for i in range(1, N + 1):
    dt = (time[i]-time[i-1]) / N_add
    c_A = c_wave(dt, con[i - 1][0], con[i - 1][1], con[i-1][2], N_add, k_vars, k_min, k_max, n)[0]
    c_B = c_wave(dt, con[i - 1][0], con[i - 1][1], con[i-1][2], N_add, k_vars, k_min, k_max, n)[1]
    c_C = c_wave(dt, con[i - 1][0], con[i - 1][1], con[i-1][2], N_add, k_vars, k_min, k_max, n)[2]
    print(i, N)
    F += (con[i][0] - c_A) ** 2 + (con[i][1] - c_B) ** 2 + (con[i][2] - c_C) ** 2



int_F = F.to_qubo()
res = qv.sim.anneal_qubo(F, num_anneals=5000)

F_solution = res.best.state

k_all = [[i[0].value(F_solution) for i in j] for j in k]

print(k_all)

k1 = [val**(1/(i+1)) for i, val in enumerate(k_all[0][1:])]
k2 = [val**(1/(i+1)) for i, val in enumerate([   asdf[0] for asdf in k_all[1:]  ])]

print('k_1=', k1)
print('k_2=', k2)

print('k_1=', sum(k1)/len(k1))
print('k_2=', sum(k2)/len(k2))

# n=nadd=3
# прошлое

# k_1= [0.59375, 1.2287252032085938, 1.5875302244850011, 1.861212084850226, 1.6437523520515278, 1.7817974539812453, 1.9622097612957417]
# k_2= [0.828125, 1.733178077982756, 1.2602285727408407, 1.1892388689867741, 1.5157174338093817, 1.6983813576847893, 1.4859942946727085]
# k_1= 1.5227110114103337
# k_2= 1.3872662294110358

# второе 

# [[1.0, 0.125, 0.015625, 2.00146484375, 4.00018310546875, 16.000015258789062, 3.814697265625e-06, 96.00000011920929], [1.0625, 1.01171875, 0.001953125, 4.00018310546875, 3.0517578125e-05, 32.00000190734863, 4.76837158203125e-07], [2.0078125, 4.0009765625, 8.0001220703125, 8.000022888183594, 48.000000953674316, 32.00000035762787], [4.0009765625, 4.00018310546875, 24.00000762939453, 32.00000190734863, 4.76837158203125e-07], [8.0001220703125, 16.000015258789062, 32.00000190734863, 32.00000035762787], [3.0517578125e-05, 16.00000286102295, 4.76837158203125e-07], [48.000000953674316, 96.00000011920929], [64.00000023841858]]
# k_1= [0.125, 0.125, 1.2602285727408407, 1.4142297464851072, 1.741101458680807, 0.12500000000000003, 1.9194712199179094]
# k_2= [1.0625, 1.4169730060943293, 1.5875302244850011, 1.6817992460012314, 0.12499999999999999, 1.9063685923065632, 1.8114473294918372]
# k_1= 0.9585758568320949
# k_2= 1.3702311997684231