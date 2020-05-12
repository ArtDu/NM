import argparse
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv


def get_points(file_name):
    with open(file_name) as f:
        x = [float(num) for num in f.readline().split()]
        y = [float(num) for num in f.readline().split()]
    return x, y


def mls(k, x, y, output=None):
    x = np.array(x)
    t = np.array([[i ** j for i in x] for j in reversed(range(k))])
    t_trans = np.transpose(t)
    g = t @ t_trans
    a = inv(g) @ t @ y
    if output is not None:
        with open(output, 'a') as f:
            coef = t @ y
            for i in range(k - 1, -1, -1):
                f.write(f'{g[i][k - 1]}a0')
                for j in range(1, k):
                    f.write(f' + {g[i][k - j - 1]}a{j}')
                f.write(f' = {coef[i]}\n')
    return a


def cal_err(x, y, coef2, coef3):
    y_err2 = []
    y_err3 = []
    for i in x:
        y_err2.append(np.polyval(coef2, i))
        y_err3.append(np.polyval(coef3, i))
    err1 = sum([(y_err2[idx] - i) ** 2 for idx, i in enumerate(y)])
    err2 = sum([(y_err3[idx] - i) ** 2 for idx, i in enumerate(y)])
    return err1, err2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input test file')
    parser.add_argument('--output', required=True, help='File for answer')
    args = parser.parse_args()

    x, y = get_points(args.input)

    with open(args.output, 'w') as f:
        pass
    coef2 = mls(2, x, y, args.output)
    coef3 = mls(3, x, y, args.output)

    # print(coef2, coef3)

    area = np.linspace(x[0], x[-1], 50)
    y_coef2 = [np.polyval(coef2, i) for i in area]
    y_coef3 = [np.polyval(coef3, i) for i in area]

    plt.scatter(x, y, color='r')
    plt.plot(area, y_coef2, color='b')
    plt.plot(area, y_coef3, color='g')
    plt.show()

    with open(args.output, 'a') as f:
        err1, err2 = cal_err(x, y, coef2, coef3)
        f.write(f'\nax + b:\n Coefficients: {list(reversed(coef2))},\n Error: {round(err1, 5)}\n')
        f.write(f'ax^2 + bx + c:\n Coefficients: {list(reversed(coef3))},\n Error {round(err2, 5)}\n')


if __name__ == "__main__":
    main()
