import numpy as np
from two_stock_fds import fds_2s


def main():
    K = 75
    sigma = .3
    r = .01
    n_ind = 2000
    nstep_price = 100
    T = 1
    nstep_time = int((sigma * nstep_price) ** 2)
    time = np.linspace(0, 1, nstep_time)
    max_s1 = .5 * K * 3
    stock = np.linspace(0, max_s1, nstep_price)
    #ret = fds_1(r, time, .5 * K, stock, sigma)
    # print(ret.shape)
    #print(np.sum(ret < 0))


if __name__ == '__main__':
    main()
