import numpy as np
import matplotlib.pyplot as plt
from runge_euler import *
from mrs import *


def mrs_real_res(x):
    return -4 * np.sin(2 * x) + 4 * x


# Zdefiniować tu funkcję y' = f(x,y(x))
def f_xy(x, y):
    return 6 * y * (np.sin(3 * x)) + 12 * np.sin(3 * x) * np.cos(3 * x)


def real_fx(x):
    return np.exp(-2 * np.cos(3 * x)) - 2 * np.cos(3 * x) + 1


# --------------------------
# Plotting functions
# --------------------------
# Function for plotting zad 1
def plot_solution(real_x, real_y, sol_x, sol_y, type):
    plt.grid(True)
    plt.plot(real_x, real_y, c="olive")
    plt.plot(sol_x, sol_y, c="red")
    plt.savefig("./" + type + "_plots/intervals=" + str(len(sol_x) - 1) + ".pdf")
    plt.close()


# Function for plotting zad 2
def plot_differential(computed_x, computed_y):
    plt.grid(True)
    real_ys = [mrs_real_res(x) for x in computed_x]
    plt.plot(computed_x, real_ys, c="olive")
    plt.plot(computed_x, computed_y, c="red")
    plt.savefig("./MRS_plots/intervals=" + str(len(computed_x) - 1) + ".pdf")
    plt.close()


def save_diffs(data, filename):
    file = open(filename, "w")
    for entry in data:
        file.write(str(entry[0]) + ';' + str(entry[1]) + '\n')


def MRS_eval():
    x0 = 0
    xk = ((2 * np.pi + 4) / 2)
    y0 = 0
    yk = mrs_real_res(xk)
    diffs = []
    for n in range(3, 50):
        diff_eq = mrs(x0, xk, n, y0, yk)
        diffs.append((n, max(abs(diff_eq[1] - yk))))
        plot_differential(diff_eq[0], diff_eq[1])

    save_diffs(diffs, "./MRS_diffs.csv")


def eul_rk_eval():
    x0 = np.pi / 6
    xk = 3 * np.pi / 2
    y0 = real_fx(x0)
    diffs_eul = []
    diffs_rk4 = []
    diffs_rk2 = []

    for n in range(3, 50):
        eul_res = euler(x0, xk, n, y0, f_xy)
        plot_solution(eul_res[0], eul_res[1], eul_res[0], real_fx(eul_res[0]), "euler")
        diffs_eul.append((n, max(abs(eul_res[1] - real_fx(eul_res[0])))))

        rk2_res = runge_kutty_2(x0, xk, n, y0, f_xy)
        plot_solution(rk2_res[0], rk2_res[1], rk2_res[0], real_fx(rk2_res[0]), "rk2")
        diffs_rk2.append((n, max(abs(rk2_res[1] - real_fx(rk2_res[0])))))

        rk4_res = runge_kutty_4(x0, xk, n, y0, f_xy)
        plot_solution(rk4_res[0], rk4_res[1], rk4_res[0], real_fx(rk4_res[0]), "rk4")
        diffs_rk4.append((n, max(abs(rk4_res[1] - real_fx(rk4_res[0])))))

    save_diffs(diffs_rk2, './rk2_diffs.csv')
    save_diffs(diffs_rk4, './rk4_diffs.csv')
    save_diffs(diffs_eul, './eul_diffs.csv')


def main():
    MRS_eval()
    eul_rk_eval()
    # print(real_fx(np.pi/6))
    # print(mrs_real_res((2*np.pi+2)/2))

if __name__ == "__main__":
    main()
