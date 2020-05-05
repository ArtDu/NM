import matplotlib.pyplot as plt
import numpy as np

def draw(x, y, name):
    plt.plot(x, y)
    plt.grid(True)
    plt.savefig(name)
    plt.cla()

def f1(x):
    return 8 / (x**2 + 4)


def f2(x):
    return ((4 - (x - 1)**2)**0.5 + 1)


def main():
    x = np.linspace(-5, 5, 10000)
    y = [f1(i) for i in x]
    y1 = [f2(i) for i in x]
    y2 = [-f2(i) + 2 for i in x]

    fig, ax = plt.subplots()
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')

    plt.plot(x, y, color='b')
    plt.plot(x, y1, color='g')
    plt.plot(x, y2, color='g')
    plt.grid(True)
    plt.savefig('p2_f')

    # для написания зависимости, лень было писать код) 
    # x = [1,2,3,4,5,6,7,8,9,10,11]
    # y = [2.461437097899386,
    #     0.5267379512665964,
    #     0.20098256531666814,
    #     0.0711066457333385,
    #     0.02596390155361563,
    #     0.009381171453750633,
    #     0.0034032981297863448,
    #     0.0012329446679830273,
    #     0.00044691033381895074,
    #     0.00016196458192323533,
    #     5.8701759597405055e-05]
    # draw(x, y, 'p2_dep_iter')

    # x = [1,2,3,4]
    # y = [0.5488721804511276,
    #     0.08230166419036067,
    #     0.0019380266142916547,
    #     1.111192574398956e-06]
    # draw(x, y, 'p2_dep_newton')

if __name__ == "__main__":
    main()