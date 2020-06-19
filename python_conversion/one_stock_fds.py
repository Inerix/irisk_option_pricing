import numpy as np


def fds_1s(r, time, strike, S, sigma):
    l_t = len(time)
    dt = time[1]
    ds = S[1]
    T = 1

    V = np.zeros((l_t, len(S)))
    V[1, :] = np.maximum(S - strike, np.zeros(len(S)))
    V[:, 1] = np.maximum(
        S[0] - strike * np.exp(-1 * r * (T - time)), np.zeros(l_t))
    V[:, (len(S) - 1)] = np.maximum(S[len(S) - 1] -
                                    strike * np.exp(-1 * r * time), np.zeros(len(time)))

    for t in range(l_t - 1):
        for i in range(1, (len(S) - 2)):
            V[t + 1, i] = (
                V[t, i] * (1 - (sigma ** 2 * S[i] ** 2 * dt) / (ds ** 2) - r * dt) +
                V[t, i + 1] * (dt / 2) * (((sigma ** 2 * S[i] ** 2) / (ds ** 2)) + (r * S[i]) / (ds)) +
                V[t, i - 1] * ((1 / 2) * ((sigma ** 2 * S[i] ** 2 * dt) / (ds ** 2)) - (r * S[i] * dt) / (2 * ds)))
    return V
