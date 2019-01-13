import numpy as np
import matplotlib.pyplot as plt


# Zdefiniować tu funkcję y' = f(x,y(x))
def f_xy(x, y):
    return 6 * y * (np.sin(3 * x)) + 12 * np.sin(3 * x) * np.cos(3 * x)


def real_fx(x):
    return np.exp(-2 * np.cos(3 * x)) - 2 * np.cos(3 * x) + 1


# Przedział [x0,xk], ilość kroków n, wartość początkowa y(x0) to y0
# Funkcja zwraca tuplę (x,y)
def euler(x0, xk, n, y0):
    y_i = y0
    ls = np.linspace(x0, xk, n)
    h = ls[1] - ls[0]
    y_s = []
    for x_i in ls:
        y_s.append(y_i)
        y_i += h * f_xy(x_i, y_i)
    return ls, y_s


# Funkcja zwraca tuplę x,y
def runge_kutty_2(x0, xk, n, y0):
    y_i = y0
    ls = np.linspace(x0, xk, n)
    h = ls[1] - ls[0]
    y_s = []
    for x_i in ls:
        y_s.append(y_i)
        k1 = h * f_xy(x_i, y_i)
        y_i += h * f_xy(x_i + h / 2, y_i + (k1 / 2))
    return ls, y_s


def runge_kutty_4(x0, xk, n, y0):
    y_i = y0
    ls = np.linspace(x0, xk, n)
    h = ls[1] - ls[0]
    y_s = []
    for x_i in ls:
        y_s.append(y_i)
        k1 = h * f_xy(x_i, y_i)
        k2 = h * f_xy(x_i + (h / 2), y_i + (k1 / 2))
        k3 = h * f_xy(x_i + (h / 2), y_i + (k2 / 2))
        k4 = h * f_xy(x_i + h, y_i + k3)
        y_i += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return ls, y_s


# --------------------------
# Różnice skończone
# Równanie to  y'' + 2xy' + 3x^2y = 3x

# Równanie to  y'' + (m ^ 2) y = (m^2)*k*x
# Równanie to  y'' + 4*y = 8*x
# ---Zad2----
# Równanie to y'' + ( m^2 * y ) = -2mcos(mx)
# Równanie to y'' + 4 * y = -4cos(2x)
# m = 2, k = 4
# p(x) funkcja przy y'
# q(x) funkcja przy y
# r(x) funkcja po prawej stronie równania (reszta)
# Warunek brzegowy -> y(0) = 0
# Warunek brzegowy 2 -> y((2pi+2/m)) = real
# --------------------------
zad = 2


def mrs_real_res(x):
    if (zad == 1):
        return -4 * np.sin(2 * x) + 4 * x
    else:
        return np.cos(2 * x) - x * np.sin(2 * x)


def p_x(x):
    if (zad == 1):
        return 0
    else:
        return 0


def q_x(x):
    if (zad == 1):
        return 4
    else:
        return 4


def r_x(x):
    if (zad == 1):
        return 16 * x
    else:
        return -4 * np.cos(2 * x)


def a_i(xi, h):
    return -2 + q_x(xi) * (h ** 2)


def b_i(xi, h):
    return 1 - ((p_x(xi) * h) / 2)


def c_i(xi, h):
    return 1 + ((p_x(xi) * h) / 2)


def mrs(x0, xk, n, u0, un):
    y_s = []
    ls = np.linspace(x0, xk, n)
    h = ls[1] - ls[0]
    z = []

    # First value
    z.append((h ** 2) * r_x(ls[0]) - b_i(ls[0], h) * u0)

    def beta(i):
        return b_i(ls[i], h) / alpha(i - 1)

    def alpha(i):
        if i == 0:
            return a_i(i, h)
        else:
            return a_i(i, h) - beta(i) * gamma(i - 1)

    def gamma(x):
        return c_i(x, h)

    for i in range(1, n - 1):
        z.append((h ** 2) * r_x(ls[i]) - beta(i) * z[i - 1])

    # Check if this is good
    # z.append(r_x(ls[-1]) - c_i(ls[-1], h) * un)
    # Maybe this
    z.append((h ** 2) * r_x(ls[-1]) - c_i(ls[-1], h) * un)

    y_s.append(z[-1] / alpha(len(z) - 1))
    for i in range(len(z) - 2, -1, -1):
        y_s.append((z[i] - gamma(i) * y_s[-1]) / alpha(i))

    return ls, y_s


# --------------------------
# Plotting functions
# --------------------------
# Function for plotting zad 1
def plot_solution(real_x, real_y, sol_x, sol_y, type):
    plt.grid(True)
    plt.plot(real_x, real_y, c="olive")
    plt.plot(sol_x, sol_y, c="red")
    plt.savefig("./" + type + "_plots/step=" + str(real_x[1] - real_x[0]) + ".pdf")
    plt.show()


# Function for plotting zad 2
def plot_differential(computed_x, computed_y):
    plt.grid(True)
    real_ys = [mrs_real_res(x) for x in computed_x]
    plt.plot(computed_x, real_ys, c="olive")
    plt.plot(computed_x, computed_y, c="red")
    plt.show()


def main():
    x0 = np.pi / 6
    xk = 3 * np.pi / 2
    n = 20  # Steps
    y0 = real_fx(x0)

    eul_res = euler(x0, xk, n, y0)
    rk2_res = runge_kutty_2(x0, xk, n, y0)
    rk4_res = runge_kutty_4(x0, xk, n, y0)

    x0 = 0

    if (zad == 1):
        xk = ((2 * np.pi + 4) / 2)
        y0 = 0
    else:
        xk = ((2 * np.pi + 2) / 2)
        y0 = 1

    n = 50  # Steps

    yk = mrs_real_res(xk)

    diff_eq = mrs(x0, xk, n, y0, yk)

    plot_solution(eul_res[0], eul_res[1], eul_res[0], real_fx(eul_res[0]), "euler")
    plot_solution(rk2_res[0], rk2_res[1], rk2_res[0], real_fx(rk2_res[0]), "rk2")
    plot_solution(rk4_res[0], rk4_res[1], rk4_res[0], real_fx(rk4_res[0]), "rk4")
    # plot_differential(diff_eq[0], diff_eq[1])


if __name__ == "__main__":
    main()
