import json
import scipy.special
import qubovert as qv
from typing import List
import matplotlib.pyplot as plt


def k_pows_list(k: qv._pubo.PUBO, n_pows: int) -> List[qv._pubo.PUBO]:
    temp = qv.PUBO()
    temp[()] = 1
    return [temp] + [k**j for j in range(1, n_pows)]


def get_mininization_func(cons: list, times: list, k_pows: list, N_exp: int, N_add: int) -> qv._pubo.PUBO:
    F = 0

    for i in range(N_exp - 1):  # !!!!!!!!
        delta = cons[i + 1][0]
        dt = (times[i + 1] - times[i]) / (N_add + 1)
        for j in range(N_add + 2):  # +1 !!!
            C = scipy.special.binom(N_add + 1, j)
            delta -= C * cons[i][0] * (dt ** j) * k_pows[j]

        F += delta ** 2

    return F


def main():
    with open('first_order.json') as f:
        data = json.load(f)

        time = data['time']
        con = data['concentrations']

        N = len(time)  # number of experimental points
        n = 10  # number of binary digits, representing k
        # number of additional dots between each pair of experimental dots
        # todo: choose with |1-kdt*dT|<1 && dt< min(dT_i)
        N_add = 4

        k_bin = [qv.PUBO.create_var('k%d' % i) for i in range(n)]

        k_dec = -sum((0.5 ** i) * k_bin[i] for i in range(n))
        k_pows = k_pows_list(k_dec, N_add+2)

        F = get_mininization_func(con, time, k_pows, N, N_add)
        res = qv.sim.anneal_pubo(F, num_anneals=1000)

        plt.hist([sample.value for sample in res])
        plt.xlabel("$F$")
        plt.ylabel("counts")
        plt.title("Values of $F$ sampled")
        plt.show()


if __name__ == "__main__":
    main()
