import json
import numpy as np
import scipy.special
import qubovert as qv

with open('Data/ChemKin/first_order.json') as f:
    data = json.load(f)

time = np.array(data['time'])
con = np.array(data['concentrations'])
N = len(time)
M = len(con[0])
K = 1  # number of k's
n = 10  # digits in k  + 1

N_add = 4# number of additional dots between each pair
          # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i) !!!  
# dT = 1
# r = int(dT/dt)
k_min = 0.125
k_max = 2

k_vars = [qv.PUBO.create_var('k%d' % i) for i in range(n)] #k = [kmin,kmax-step]
k = -sum((0.5 ** i) * k_vars[i] for i in range(n))
print(k_vars)
print("type(k) = ", type(k))
#for j in range(3):
#    print(k**j)      #doesn't work :(

temp = qv.PUBO()
temp[()] = 1
k_pows = []
for i in range(N_add+2):
    #print("i=",i,"k^i=",temp)
    print("============temp============")
    print(temp)
    print("============================")
    k_pows.append(temp)
    temp = temp*k

print("---------------k_pows-------------")
print("k_pows[0]: ",k_pows[0])
print("k_pows[1]: ",k_pows[1])
print("k_pows[2]: ",k_pows[2])
print("----------------------------------")


F = 0

for i in range(N-1): #!!!!!!!!
    delta = con[i+1][0]
    dt = (time[i+1]-time[i]) / (N_add + 1)
    #print("dt =",dt)
    for j in range(N_add+2): # +1 !!!
        C = scipy.special.binom(N_add+1, j)
        
        #delta -= C * con[i - 1][0] * (dt * k)**j
        delta -= C * con[i][0] * (dt**j) * k_pows[j]
    #print(i, delta)
    F += delta ** 2

print("F =")
print(F)

print()
int_F = F.to_pubo()
print("F vars:", F.variables)
print("int_F vars:", int_F.variables)

"""
int_F_solution = int_F.solve_bruteforce()
print("bruteforce solution:", int_F_solution)
F_solution = F.convert_solution(int_F_solution)
print(F_solution)
print()
print("k =", k.value(F_solution))
print("F(k) =", F.value(F_solution))
"""
print("degree =",F.degree) 

print()
print("anneal_pubo solution:")
res = qv.sim.anneal_pubo(F, num_anneals=1000)
print(res.best)
print()

F_solution = res.best.state
print("F_solution =", F_solution)
print("k =", k.value(F_solution))
print("F(k) =", F.value(F_solution))



import matplotlib.pyplot as plt
plt.hist([sample.value for sample in res])
plt.xlabel("$F$")
plt.ylabel("counts")
plt.title("Values of $F$ sampled")
plt.show()
