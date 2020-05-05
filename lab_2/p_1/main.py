import math
import numpy as np
import matplotlib.pyplot as plt


# 2^x - x ^ 2 - 0.5 = 0
def f(x):
    return 2 ** x - x ** 2 - 0.5

def df(x):
    return 2 ** x * math.log(2) - 2 * x

def ddf(x):
    return 2 ** x * math.log(2) ** 2 - 2

def phi(x):
    return math.sqrt(2 ** x - 0.5)


def dphi(x):
    return (2 ** x * math.log(2)) / (2 * math.sqrt(2 ** x - 0.5))


def simpleIteration(phi, dphi, a, b, eps=0.0001):
    q = max(abs(dphi(a)), abs(dphi(b)))
    x = (a + b) / 2
    k = 0
    go = True

    while go:
        k += 1
        x_cur = phi(x)
        # ,q/(1-q)*|x_cur - x|: {(q * abs(x_cur - x) / (1 - q))}
        print(f'x: {round(x_cur, 4)}, k: {k}')
        if (q * abs(x_cur - x) / (1 - q)) <= eps:
            go = False

        x = x_cur
        if k == 100:
            break


def newton(f, df, x0, eps=0.0001):
    x = x0
    k = 0
    go = True
    while go:
        k += 1
        x_cur = x - f(x) / df(x)
        # , |x_cur - x|: {abs(x_cur - x)}
        print(f'x: {round(x_cur, 4)}, k: {k}')
        if abs(x_cur - x) <= eps:
            go = False

        x = x_cur


def show(f, df, x, file = None, step = 0.5, ddf = None, label = None):
    X = np.arange(x[0], x[-1], step)
    Y = [f(i) for i in X]
    dY = [df(i) for i in X]

    if ddf:
        ddY = [ddf(i) for i in X]

    fig, axis = plt.subplots()
    axis.plot(X, Y, label=label[0])
    axis.plot(X, dY, label=label[1])

    if ddf:
        axis.plot(X, ddY, label=label[2])

    axis.legend(loc='upper right')
    axis.grid()

    if file:
        fig.savefig(file)
        print(f'File {file} was saved correctly')
        plt.close(fig)

    plt.show()


if __name__ == '__main__':
    print("My function is: 2^x - x ^ 2 - 0.5")
    eps = float(input("Enter epsilon:\n"));

    print("Simple iteration:")
    simpleIteration(phi, dphi, 1, 2, eps)
    print("Newton method:")
    newton(f, df, 1.5, eps)
    show(f, df, [0, 2], step=0.1, ddf=ddf, label=['f','df','dff'])
    show(phi, dphi, [0, 2], step=0.1, label=['phi', 'dphi'])
