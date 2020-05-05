import matplotlib.pyplot as plt
import numpy as np
import Lw2_1


def draw(x, y, name):
    plt.plot(x, y)
    plt.grid(True)
    plt.savefig(name)
    plt.cla()


def main():
    a, b = (1, 3)

    x = np.linspace(-5, 5, 1000)
    y = [Lw2_1.Solver.f(i) for i in x]
    draw(x, y, 'p1_f')

    x = [1,2,3,4,5,6]
    y = [0.4376792356304313,
        0.09868531531864919,
        0.013477821811890169,
        0.0016260850279907106,
        0.000192826088396256,
        2.281817613340577e-05]
    draw(x, y, 'p1_dep_iter')

    x = [1,2,3,4]
    y = [0.37104878347507464,
        0.04667937242463971,
        0.000540410722716711,
        7.66827679132831e-08]
    draw(x, y, 'p1_dep_newton')


if __name__ == "__main__":
    main()