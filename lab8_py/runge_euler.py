import numpy as np


# Przedział [x0,xk], ilość kroków n, wartość początkowa y(x0) to y0
# Funkcja zwraca tuplę (x,y)
def euler(x0, xk, n, y0, f_xy):
    y_i = y0
    ls = np.linspace(x0, xk, n)
    h = ls[1] - ls[0]
    y_s = []
    for x_i in ls:
        y_s.append(y_i)
        y_i += h * f_xy(x_i, y_i)
    return ls, y_s


# Funkcja zwraca tuplę x,y
def runge_kutty_2(x0, xk, n, y0, f_xy):
    y_i = y0
    ls = np.linspace(x0, xk, n)
    h = ls[1] - ls[0]
    y_s = []
    for x_i in ls:
        y_s.append(y_i)
        k1 = h * f_xy(x_i, y_i)
        y_i += h * f_xy(x_i + h / 2, y_i + (k1 / 2))
    return ls, y_s


def runge_kutty_4(x0, xk, n, y0, f_xy):
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
