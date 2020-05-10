import argparse
import numpy as np


def lagrange(points):
    def inner(x):
        cnt_points = len(points)
        l = []
        omega_apostrophe = []
        for i in range(cnt_points):
            l.append(np.prod([(x - points[j][0]) / (points[i][0] - points[j][0])
                              for j in range(cnt_points) if i != j]))
            omega_apostrophe.append(round(np.prod([(points[i][0] - points[j][0])
                                                   for j in range(cnt_points) if i != j]), 3))
        result = sum([points[idx][1] * i for idx, i in enumerate(l)])
        f_div_omega_ap = [round(points[idx][1] / i, 5) for idx, i in enumerate(omega_apostrophe)]
        x_minus_xi = [round(x - xi_fi[0], 1) for xi_fi in points]
        table_i = [i for i in range(cnt_points)]
        table_xi = [xi_fi[0] for xi_fi in points]
        table_fi = [round(xi_fi[1], 5) for xi_fi in points]
        return result, table_i, table_xi, table_fi, omega_apostrophe, f_div_omega_ap, x_minus_xi

    return inner


def newton(points):
    def inner(x):
        cnt_points = len(points)
        x_point, y_point = zip(*points)
        coef = [y_point[0]]
        for j, shift in zip(reversed(range(1, cnt_points)), range(cnt_points)):
            tmp = []
            for l, r in zip(range(j), range(1, j + 1)):
                num1 = y_point[l] - y_point[r]
                num2 = x_point[l] - x_point[r + shift]
                tmp.append(num1 / num2)
            y_point = tmp
            coef.append(y_point[0])

        res = 0
        for i in range(cnt_points):
            res += coef[i] * np.prod([x - x_point[j] for j in range(i)])
        P = f'{coef[0]} + {round(coef[1],5)}(x - {x_point[0]}) + ' \
            f'{round(coef[2],5)}(x - {x_point[0]})(x - {x_point[1]}) + ' \
            f'{round(coef[3],5)}(x - {x_point[0]})(x - {x_point[1]})(x - {x_point[2]}).\n'

        return res, P

    return inner


def get_points(file_name):
    with open(file_name) as f:
        x = [float(num) * np.pi for num in f.readline().split()]
        # x = [float(num) for num in f.readline().split()]
        points = [(round(i,5), round(np.sin(i), 5)) for i in x]
        x_test = float(f.readline()) * np.pi
        # x_test = float(f.readline())
    return points, x_test


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input test file')
    parser.add_argument('--output', required=True, help='File for answer')
    args = parser.parse_args()

    points, x_test = get_points(args.input)
    lg = lagrange(points)
    nt = newton(points)

    with open(args.output, 'w') as f:
        tmp = round(x_test, 3)
        res1, table_i, table_xi, table_fi, omega_apostrophe, f_div_omega_ap, x_minus_xi = lg(tmp)
        true_val = np.sin(tmp)
        # true_val = np.log(tmp)
        eps1 = abs(res1 - true_val)
        f.write(f'Lagrange:\n L({tmp}) = {res1}\n')
        f.write(f' y({tmp}): {round(true_val, 5)}\n')
        f.write(f' Eps = {round(eps1, 5)}\n')
        f.write(f' i: {table_i}\n')
        f.write(f' xi: {table_xi}\n')
        f.write(f' fi: {table_fi}\n')
        f.write(f' omega\'(xi): {omega_apostrophe}\n')
        f.write(f' fi/omega\'(xi): {f_div_omega_ap}\n')
        f.write(f' X* - xi: {x_minus_xi}\n')

        res2, P = nt(tmp)
        eps2 = abs(res2 - true_val)
        f.write(f'Newton:\n N({tmp}) = {res2}\n ')
        f.write(P)
        f.write(f' Eps = {round(eps2, 5)}\n')
        f.write(f'Orig:\n Sin({tmp}) = {true_val}\n')


if __name__ == "__main__":
    main()
