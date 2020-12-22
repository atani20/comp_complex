import numpy as np
from matplotlib import pyplot as plt
import os
from kaucherpy_master.src.kaucher.Kaucher import *


def sti_vector(x_inf, x_sup):
    return np.hstack([-x_inf, x_sup])


def back_sti_vector(sti_vec):
    half = sti_vec.shape[0] // 2
    x_inf = -sti_vec[:half]
    x_sup = sti_vec[half:]
    return (x_inf, x_sup)


def sti_dot_mtx(a_mtx):
    a_plus = np.zeros(a_mtx.shape)
    a_minus = np.zeros(a_mtx.shape)
    for i in range(a_mtx.shape[0]):
        for j in range(a_mtx.shape[1]):
            if (a_mtx[i][j] > 0):
                a_plus[i][j] = a_mtx[i][j]
            else:
                a_minus[i][j] = -a_mtx[i][j]
    return np.block([[a_plus, a_minus], [a_minus, a_plus]])


def compute_G(C_mtx, cur_point, d_sti):
    intervaled_inf, intervaled_sup = back_sti_vector(cur_point)
    y_interval = [Kaucher(intervaled_inf[i], intervaled_sup[i]) for i in range(intervaled_inf.shape[0])]
    dot_prod = [sum([C_mtx[i][j] * y_interval[j] for j in range(len(y_interval))]) for i in range(len(C_mtx))]
    prod_lower = np.array([comp.lower for comp in dot_prod])
    prod_upper = np.array([comp.upper for comp in dot_prod])
    return sti_vector(prod_lower, prod_upper) - d_sti


def compute_x0(mid_C, sti_d):
    return np.linalg.solve(mid_C, sti_d)


def compute_D(D, i, j, a_ij, b_inf, b_sup):
    n = D.shape[0] // 2;

    ainf = a_ij.lower
    asup = a_ij.upper

    # determine case (1 of 16)
    k = 0
    m = 0
    if ainf * asup > 0:
        k = 0 if ainf > 0 else 2
    else:
        k = 1 if ainf < asup else 3

    if b_inf * b_sup > 0:
        m = 1 if b_inf > 0 else 3
    else:
        m = 2 if b_inf <= b_sup else 4

    case = 4 * k + m
    if case == 1:
        D[i, j] = ainf
        D[i + n, j + n] = asup
    elif case == 2:
        D[i, j] = asup
        D[i + n, j + n] = asup
    elif case == 3:
        D[i, j] = asup
        D[i + n, j + n] = ainf
    elif case == 4:
        D[i, j] = ainf
        D[i + n, j + n] = ainf
    elif case == 5:
        D[i, j + n] = ainf
        D[i + n, j + n] = asup
    elif case == 6:
        if ainf * b_sup < asup * b_inf:
            D[i, j + n] = ainf
        else:
            D[i, j] = asup
        if ainf * b_inf > asup * b_sup:
            D[i + n, j] = ainf
        else:
            D[i + n, j + n] = asup
    elif case == 7:
        D[i, j] = asup
        D[i + n, j] = ainf
    elif case == 8:
        pass  # zeroes are already filled
    elif case == 9:
        D[i, j + n] = ainf
        D[i + n, j] = asup
    elif case == 10:
        D[i, j + n] = ainf
        D[i + n, j] = ainf
    elif case == 11:
        D[i, j + n] = asup
        D[i + n, j] = ainf
    elif case == 12:
        D[i, j + n] = asup
        D[i + n, j] = asup
    elif case == 13:
        D[i, j] = ainf
        D[i + n, j] = asup
    elif case == 14:
        pass  # zeros are already filled
    elif case == 15:
        D[i, j + n] = asup
        D[i + n, j + n] = ainf
    elif case == 16:
        if ainf * b_inf > asup * b_sup:
            D[i, j] = ainf
        else:
            D[i, j + n] = -asup
        if ainf * b_sup < asup * b_inf:
            D[i + n, j + n] = ainf
        else:
            D[i + n, j] = asup
    return D


def subdiff2(A_inf, A_sup, inf_b, sup_b):
    # params
    learning_rate = 0.8
    accuracy = 1e-8
    max_iter = 40

    n = A_inf.shape[0]  # dim

    # compute x0
    A_mid = np.array([[(A_inf[i, j] + A_sup[i, j]) / 2 for j in range(n)] for i in range(n)])
    A = [[Kaucher(A_inf[i, j], A_sup[i, j]) for j in range(n)] for i in range(n)]
    C_sti = sti_dot_mtx(A_mid)  # for x0
    d_sti = sti_vector(inf_b, sup_b)  # for x0
    # cur_x = compute_x0(C_sti, d_sti)  # x0: (mid C)~ x0 = sti(b)
    cur_x = np.zeros(d_sti.shape[0])  # x0 = 0

    prev_x = cur_x
    started = False  # to enter the first iteration

    # do until method converges
    cur_iter = 0
    worklist = [cur_x]  # return each step
    while not started or np.linalg.norm(cur_x - prev_x) > accuracy:
        started = True
        cur_iter += 1
        if (cur_iter > max_iter):
            print("Too many iterations")
            break
        prev_x = cur_x  # update after the previous step
        D = np.zeros((2 * n, 2 * n))  # subgrad mtx
        for i in range(n):
            for j in range(n):
                g = A[i][j]
                h_inf = -prev_x[j]
                h_sup = prev_x[j + n]
                D = compute_D(D, i, j, g, h_inf, h_sup)
        G_part = compute_G(A, prev_x, d_sti)
        dx = np.linalg.solve(D, -G_part)
        cur_x = prev_x + learning_rate * dx
        worklist.append(cur_x)
    return back_sti_vector(cur_x), worklist, cur_iter

a_slae_inf = np.array([[1, 1], [0, 0.1]])
a_slae_sup = np.array([[1, 1], [0, 0.2]])

inf_b = np.array([4.1, 0.1])
sup_b = np.array([3.9, 0.4])


(x_inf, x_sup), x_seq, iter_count = subdiff2(a_slae_inf, a_slae_sup, inf_b, sup_b)
print(f"Iters: {iter_count}")
print(f"X inf: {x_inf}")
print(f"X sup: {x_sup}")

