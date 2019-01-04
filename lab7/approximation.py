import numpy as np
import matplotlib.pyplot as plt
import os
import math


def f(x):
    return np.exp(-3 * np.sin(2 * x))


def get_d(xpts, degree):
    n = len(xpts)
    D = np.zeros((n, degree))
    for i in range(n):
        for j in range(degree):
            D[i][j] = xpts[i] ** j

    return D


def get_coeffs(nodes, degree):
    xs = [node[0] for node in nodes]
    ys = [node[1] for node in nodes]
    D = get_d(xs, degree)
    Dt = np.transpose(D)
    DtD = np.dot(Dt, D)
    invDtD = np.linalg.inv(DtD)  # we inverse the matrix
    return np.dot(invDtD, (np.dot(Dt, ys)))  # where the magic happens


# Returns a tuple with (ak,bk)
def get_trig_coeff(nodes, k):
    akl = []
    bkl = []
    m = len(nodes)
    for i in range(m):
        akl.append(nodes[i][1] * np.cos(k * nodes[i][0]))
        bkl.append(nodes[i][1] * np.sin(k * nodes[i][0]))
    return (2 * sum(akl) / m, 2 * sum(bkl) / m)


def approximate_trig(x, n, nodes):  # N is for degree of the approximation
    a0 = get_trig_coeff(nodes, 0)[0]
    an = get_trig_coeff(nodes, n)[0]
    res = 0
    for k in range(1, n):
        ak, bk = get_trig_coeff(nodes, k)
        res += ak * np.cos(k * x)
        res += bk * np.sin(k * x)
    res += a0 / 2
    res += an * np.cos(n * x)
    return res


def scatter_nodes(nodes, cheby_nodes, result):
    nxs = [p[0] for p in nodes]
    nys = [p[1] for p in nodes]
    nchxs = [p[0] for p in cheby_nodes]
    nchys = [p[1] for p in cheby_nodes]
    plotxs = [p[0] for p in result]
    plotys = [p[1] for p in result]

    plt.plot(plotxs, plotys, c="olive")
    plt.scatter(nxs, nys, c="black")
    plt.scatter(nchxs, nchys, c="orange")
    plt.show()


def plot_result(original, nodes, result, cheby_nodes, cheby_result, degree):
    origxs = [p[0] for p in original]
    origys = [p[1] for p in original]
    nodesx = [p[0] for p in nodes]
    nodesy = [p[1] for p in nodes]
    # resx = [p[0] for p in result]
    resy = [p[1] for p in result]
    chresy = [p[1] for p in cheby_result]
    chnodesy = [p[1] for p in cheby_nodes]
    chnodesx = [p[0] for p in cheby_nodes]

    plt.plot(origxs, origys, c='olive')
    plt.plot(origxs, resy, c='red')
    plt.scatter(nodesx, nodesy, c='black')
    plt.grid(True)
    dir = "./plots"
    if not os.path.exists(dir):
        os.mkdir(dir)
    dir += '/degree=' + str(degree) + '/'
    if not os.path.exists(dir):
        os.mkdir(dir)
    plt.savefig(dir + "EVEN_nodes=" + str(len(nodes)) + ".pdf")
    plt.close()
    # Now cheby part
    plt.grid(True)
    plt.plot(origxs, origys, c='olive')
    plt.plot(origxs, chresy, c='blue')
    plt.scatter(chnodesx, chnodesy, c='black')
    dir = "./plots"
    if not os.path.exists(dir):
        os.mkdir(dir)
    dir += '/degree=' + str(degree) + '/'
    if not os.path.exists(dir):
        os.mkdir(dir)
    plt.savefig(dir + "CHEBY_nodes=" + str(len(nodes)) + ".pdf")
    plt.close()


def plot_trig(original, nodes, cheby_nodes, result_trig_even, result_trig_cheby, deg):
    origxs = [p[0] for p in original]
    origys = [p[1] for p in original]
    nodesx = [p[0] for p in nodes]
    nodesy = [p[1] for p in nodes]

    evx = [p[0] for p in result_trig_even]
    evy = [p[1] for p in result_trig_even]
    # Even part plot save
    plt.plot(origxs, origys, c='olive')
    plt.plot(evx, evy, c='red')
    plt.scatter(nodesx, nodesy, c='black')
    plt.grid(True)
    dir = "./plots/trig/deg="+str(deg)+'/'
    if not os.path.exists(dir):
        os.mkdir(dir)

    plt.savefig(dir + "even_nodes=" + str(len(nodes)) + ".pdf")
    plt.close()

    # Now cheby part
    chx = [p[0] for p in result_trig_cheby]
    chy = [p[1] for p in result_trig_cheby]

    plt.grid(True)
    plt.plot(origxs, origys, c='olive')
    plt.plot(chx, chy, c='blue')
    plt.scatter([p[0] for p in cheby_nodes], [p[1] for p in cheby_nodes], c='black')

    plt.savefig(dir + "cheby_nodes=" + str(len(nodes)) + ".pdf")
    plt.close()


def approximate(coeffs, x):
    res = 0.0
    for i in range(len(coeffs)):
        res += coeffs[i] * (pow(x, i))
    return res


def cheby_x(x0, x1, n):
    result = []
    for i in range(1, n + 1, 1):
        result.append(1 / 2 * (x0 + x1) + 1 / 2 * (x1 - x0) * np.cos((2 * i - 1) * np.pi / (2 * n)))
    return result


def main():
    n = 15  # how many points we want to map
    resolution = 1000  # how big the plot's resolution will be
    degree = 10  # degree of the interpolation polynomial - 1. it means 2 is a line function 3 is quadratic etc...

    xdraw = np.linspace(0, 3 * np.pi, resolution)
    ydraw = [f(xpt) for xpt in xdraw]

    drawn = list(zip(xdraw, ydraw))

    # coeffs = get_coeffs(nodes, degree) #we get the coeffitients of a polynomial of wanted degree
    # yres = [approximate(coeffs, xpt) for xpt in xdraw]
    # result = list(zip(xdraw, yres))

    print("deg;nodes;eucl_dist")
    for n in range(4, 40, 5):  # we change the number of nodes, from 4 to 15 (n)
        for i in range(2, 10, 2):  # we change the degree of polynomial (or number of (sines + cosines / 2 ) )
            diffs = []
            cheby_diffs = []
            xnodes = np.linspace(0, 3 * np.pi, n)
            ynodes = [f(xnode) for xnode in xnodes]
            x_cheby_nodes = cheby_x(0, 3 * np.pi, n)
            y_cheby_nodes = [f(xnode) for xnode in x_cheby_nodes]

            cheby_nodes = list(zip(x_cheby_nodes, y_cheby_nodes))
            nodes = list(zip(xnodes, ynodes))

            # coeffs = get_coeffs(nodes, i)
            # coeffs_cheby = get_coeffs(cheby_nodes, i)
            #
            # yres = [approximate(coeffs, xpt) for xpt in xdraw]
            # y_cheby_res = [approximate(coeffs_cheby, xpt) for xpt in xdraw]

            y_trig_even_res = [approximate_trig(xpt, i, nodes) for xpt in xdraw]
            y_trig_cheby_res = [approximate_trig(xpt, i, cheby_nodes) for xpt in xdraw]

            # result = list(zip(xdraw, yres))
            # result_cheby = list(zip(xdraw, y_cheby_res))
            result_trig_even = list(zip(xdraw, y_trig_even_res))
            result_trig_cheby = list(zip(xdraw, y_trig_cheby_res))

            f_even = open("./EVEN_diffs.csv", "a")
            f_cheby = open("./CHEBY_diffs.csv", "a")

            # for ind in range(len(xdraw)):
            #     diffs.append((ydraw[ind] - yres[ind]) ** 2)
            # f_even.write((str(i) + ";" + str(n) + ";" + str((sum(diffs) / len(diffs)) ** (1 / 2))) + "\n")
            #
            # for ind in range(len(xdraw)):
            #     cheby_diffs.append((ydraw[ind] - y_cheby_res[ind]) ** 2)
            # f_cheby.write((str(i) + ";" + str(n) + ";" + str((sum(cheby_diffs) / len(cheby_diffs)) ** (1 / 2))) + "\n")
            #
            # plot_result(drawn, nodes, result, cheby_nodes, result_cheby, i)
            plot_trig(drawn, nodes, cheby_nodes, result_trig_even, result_trig_cheby, i)
            # scatter_nodes(nodes, cheby_nodes, drawn)


if __name__ == "__main__":
    main()
    exit(0)
