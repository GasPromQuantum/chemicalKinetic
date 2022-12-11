import json
import scipy.special
import qubovert as qv
from typing import List
import matplotlib.pyplot as plt
from utils import get_mininization_func


def k_pows_list(k: qv._pubo.PUBO, n_pows: int) -> List[qv._pubo.PUBO]:
    temp = qv.PUBO()
    temp[()] = 1
    return [temp] + [k**j for j in range(1, n_pows)]


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
        res = qv.sim.anneal_pubo(F, num_anneals=10)

        F_solution = res.best.state
        print("\nF_solution =", F_solution)
        print("k =", k_dec.value(F_solution))
        print("F(k) =", F.value(F_solution))

        plt.hist([sample.value for sample in res])
        plt.xlabel("$F$")
        plt.ylabel("counts")
        plt.title("Values of $F$ sampled")
        plt.show()

        F_qubo = F.to_qubo()
        F_qubo_solution = qv.sim.anneal_qubo(F_qubo, num_anneals=1000)
        res_qubo = F_qubo_solution.best.state
        print(type(res_qubo))
        res_pubo_from_qubo = F.convert_solution(res_qubo)

        print("\nF_solution =", res_pubo_from_qubo)
        print("k =", k_dec.value(res_pubo_from_qubo))
        print("F(k) =", F.value(res_pubo_from_qubo))

        plt.hist([sample.value for sample in F_qubo_solution])
        plt.xlabel("$F$")
        plt.ylabel("counts")
        plt.title("Values of $F$ sampled")
        plt.show()


if __name__ == "__main__":
    main()
