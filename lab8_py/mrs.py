import numpy as np
# # --------------------------
# # Różnice skończone
# # Równanie to  y'' + 2xy' + 3x^2y = 3x
#
# # Równanie to  y'' + (m ^ 2) y = (m^2)*k*x
# # Równanie to  y'' + 4*y = 8*x
# # m = 2, k = 4
# # p(x) funkcja przy y'
# # q(x) funkcja przy y
# # r(x) funkcja po prawej stronie równania (reszta)
# # Warunek brzegowy -> y(0) = 0
# # Warunek brzegowy 2 -> y((2pi+2/m)) = real
# # --------------------------
#

def p_x(x):
    return 0


def q_x(x):
    return 4


def r_x(x):
    return 16 * x


def a_i(xi, h):
    return -2 + q_x(xi) * (h ** 2)


def b_i(xi, h):
    return 1 - ((p_x(xi) * h) / 2)


def c_i(xi, h):
    return 1 + ((p_x(xi) * h) / 2)



def mrs(x0, xk, n, u0, un):
    xs = np.linspace(x0, xk, n, True)
    A = np.zeros((n, n))
    dx = xs[1] - xs[0]

    A[0][0] = 1

    for i in range(1, n-1):
        A[i][i - 1] = b_i(xs[i], dx)
        A[i][i] = a_i(xs[i], dx)
        A[i][i + 1] = c_i(xs[i], dx)

    A[n-1][n-1] = 1

    f = []

    f.append(u0)

    for i in range(1, n-1):
        f.append((dx ** 2) * r_x(xs[i]))

    f.append(un)

    # A*u = f
    # u = szukane

    u = np.linalg.solve(A,f)

    return xs, u
