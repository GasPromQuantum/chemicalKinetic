import qubovert as qv
import scipy.special

def get_mininization_func(cons: list, times: list, k_pows: list, N_exp: int, N_add: int) -> qv._pubo.PUBO:
    F = 0

    for i in range(N_exp - 1):
        delta = cons[i + 1][0]
        dt = (times[i + 1] - times[i]) / (N_add + 1)
        for j in range(N_add + 2):
            C = scipy.special.binom(N_add + 1, j)
            delta -= C * cons[i][0] * (dt ** j) * k_pows[j]

        F += delta ** 2

    return F