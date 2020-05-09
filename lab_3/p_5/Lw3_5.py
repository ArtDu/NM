import argparse
import numpy as np


# x0, xk, h1, h2
def get_points(file_name):
    with open(file_name) as f:
        a, b, h1, h2 = [float(num) for num in f.readline().split()]
    return a, b, h1, h2

def f(x):
    return x / (2 * x + 5)


def simpson(a, b, h, n):
    s = f(a) + f(b)
    for i in range(1, n, 2):
        s += 4 * f(a + i * h)
    for i in range(2, n-1, 2):
        s += 2 * f(a + i * h)
    return s * h / 3


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input test file')
    parser.add_argument('--output', required=True, help='File for answer')
    args = parser.parse_args()

    a, b, h1, h2 = get_points(args.input)

    x1 = np.linspace(a, b, int((b - a) / h1 + 1))
    x2 = np.linspace(a, b, int((b - a) / h2 + 1))

    y_trap1 = [f(i) for i in x1]
    y_trap2 = [f(i) for i in x2]

    rect1 = h1 * sum([f((i + j) / 2) for i, j in zip(x1, x1[1:])])
    rect2 = h2 * sum([f((i + j) / 2) for i, j in zip(x2, x2[1:])])

    trap1 = h1/2 * sum([i + j for i, j in zip(y_trap1[1:], y_trap1)])
    trap2 = h2/2 * sum([i + j for i, j in zip(y_trap2[1:], y_trap2)])

    simps1 = simpson(a, b, h1, int((b - a) / h1))
    simps2 = simpson(a, b, h2, int((b - a) / h2))

    with open(args.output, 'w') as fi:
        fi.write(f'Method of rectangles = {rect1}\tStep = {h1}\n')
        fi.write(f'Method of rectangles = {rect2}\tStep = {h2}\n')
        fi.write(f'Error: {round(abs(rect1 - rect2) / 3, 5)}\n')
        fi.write('='*10)
        fi.write('\n')
        fi.write(f'Method of trapeziums = {trap1}\tStep = {h1}\n')
        fi.write(f'Method of trapeziums = {trap2}\tStep = {h2}\n')
        fi.write(f'Error: {round(abs(trap1 - trap2) / 3, 5)}\n')
        fi.write('='*10)
        fi.write('\n')
        fi.write(f"Simpson's method  = {simps1}\tStep = {h1}\n")
        fi.write(f"Simpson's method  = {simps2}\tStep = {h2}\n")
        fi.write(f'Error: {round(abs(simps1 - simps2) / 15, 5)}\n')

if __name__ == "__main__":
    main()