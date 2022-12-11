import json
import numpy as np

with open('Data/ChemKin/second_order.json') as f:
    data = json.load(f)

time = np.array(data['time'])
con = np.array(data['concentrations'])
N = len(time) - 1
M = len(con[0])
K = 1
dt = 0.5

RHS = np.zeros(shape=(N * M, K))

for i in range(M):
    for j in range(N):
        RHS[i + j * M] = (con[j + 1][i] - con[j][i]) / dt


LHS = np.zeros(shape=(N * M, K))
for i in range(M):
    for j in range(N):
        LHS[i + j * M][0] = - con[j][0] * con[j][1]
LHS_t = np.transpose(LHS)

print(RHS)
print(LHS)
print(LHS_t @ RHS)
print(LHS_t @ LHS)

A = (LHS_t @ RHS)[0][0]
B = (LHS_t @ LHS)[0][0]

print(A / B)



