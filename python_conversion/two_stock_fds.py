import numpy as np
import pandas as pd

from one_stock_fds import fds_1s


def fds_2s(K, sigma=[.3, .3], cor=.3, r=.01, weights=[.5, .5], n_ind=2000, nstep_price=100):
    s1 = sigma[0]
    s2 = sigma[1]
    w1 = weights[0]
    w2 = weights[1]

    # for 1 - stock case, convergence is nstep_time >= (sigma ** 2 * I ** 2)
    # Condition can be found on Wilmott pg. 1215
    T = 1
    nstep_time = (max(sigma) * nstep_price) ** 2  # 900
    delta_t = T/nstep_time

    # Time vector increments
    time = np.linspace(0, T, delta_t)

    V = np.zeros((nstep_time + 1, nstep_price + 1, nstep_price + 1))
    # dimensions are [time, s1, s2]

    # 3 is an arbitrary number determined by Wilmott pg. 1202
    max_b1 = w1 * K * 3
    b_1 = np.linspace(0, max_b1, nstep_price + 1)

    max_b2 = w2 * K * 3
    b_2 = np.linspace(0, max_b2, nstep_price + 1)

    ds1 = b_1[1]
    ds2 = b_2[1]

    for i in range(0, len(b_1)):
        for j in range(0, len(b_2)):
            V[1, i, j] = max(b_1[i] * w1 + b_2[j] * w2 - K, 0)

    # set 4 boundary conditions

    # sets min S1
    # both of these are 0
    V[:, 1, :] = fds_1s(r, time, K - w1 * b_1[1], w2 * b_2, s2)

    # sets min S2
    V[:, :, 1] = fds_1s(r, time, K - w2 * b_2[1], w1 * b_1, s1)

    # sets max S1
    V[:, nstep_price + 1,
        :] = fds_1s(r, time, K - w1 * b_1[nstep_price + 1], w2 * b_2, s2)

    # sets max S2
    V[:, :, nstep_price +
        1] = fds_1s(r, time, K - w2 * b_2[nstep_price + 1], w1 * b_1, s1)

    # iterate through at each slice

    for t in range(nstep_time):
        for i in range(1, nstep_price):
            s1_t = b_1[i]
            s1_thing = (delta_t * (s1 * w1 * s1_t) ** 2)/(ds1 ** 2)
            r_thing_1 = (r * delta_t * s1_t) / (2 * ds1)
            for j in range(1, nstep_price):
                s2_t = b_2[j]
                s2_thing = (delta_t * (s2 * w2 * s2_t) ** 2) / ds2 ^ 2
                r_thing_2 = (r * delta_t * s2_t) / (2 * ds2)
                huge_thing = ((s1 * s2 * cor * w1 * w2 * s1_t * s2_t * delta_t) /
                              (8 * ds1 * ds2))
                V[t + 1, i, j] = (
                    V[t, i, j] * (1 - s1_thing - s2_thing - r * delta_t) +
                    V[t, i + 1, j] * (r_thing_1 + s1_thing / 2) +
                    V[t, i - 1, j] * (s1_thing / 2 - r_thing_1) +
                    V[t, i, j + 1] * (r_thing_2 + s2_thing/2) +
                    V[t, i, j - 1] * (s2_thing / 2 - r_thing_2) +
                    V[t, i + 1, j + 1] * huge_thing -
                    V[t, i - 1, j + 1] * huge_thing -
                    V[t, i + 1, j - 1] * huge_thing +
                    V[t, i - 1, j - 1] * huge_thing)
    return V[len(time) - 1, :, :]
