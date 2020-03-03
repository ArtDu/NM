import sys
import numpy as np
from math import pi, atan, cos, sin, sqrt

epsilon = 0.01


class MatrixError(Exception):
    pass


class MethodError(Exception):
    pass


def matrix_product(A, B):
    n1 = len(A)
    m1 = len(A[0])
    n2 = len(B)
    m2 = len(B[0])
    if m1 != n2:
        raise MatrixError("Matrix product error")
    result = [[0 for _ in range(m2)] for _ in range(n1)]
    for i in range(n1):
        for k in range(m2):
            for j in range(m1):
                result[i][k] += A[i][j] * B[j][k]
    return result


def QR_decompostion(A):
    n = len(A)
    v = [0 for _ in range(n)]
    A_tmp = A.copy()

    Q = [[0 if i != j else 1 for i in range(n)] for j in range(n)]

    for i in range(n - 1):
        norma = 0
        for j in range(i, n):
            norma += (A_tmp[j][i]) ** 2
        norma = sqrt(norma)
        v = [0 for _ in range(n)]
        v[i] = A_tmp[i][i] + np.sign(A_tmp[i][i]) * norma
        for j in range(i + 1, n):
            v[j] = A_tmp[j][i]

        H = [[0 if i != j else 1 for i in range(n)] for j in range(n)]

        lower = 0
        for i in range(n):
            lower += v[i] * v[i]

        for i in range(n):
            for j in range(n):
                upper = v[i] * v[j]
                H[i][j] -= 2 * (upper / lower)

        Q = matrix_product(Q, H)

        A_tmp = matrix_product(H, A_tmp)
        A_tmp = [[round(A_tmp[i][j], 2) for j in range(n)] for i in range(n)]

    return Q, A_tmp

def QR_method(A):
    A_tmp = A.copy()
    while True:

        Q, R = QR_decompostion(A_tmp)
        A_tmp = matrix_product(R, Q)

        # epsilon_real() and epsilon_complex()

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("use {} <matrix_file> <b_file>")
        exit(0)

    matrix_file = sys.argv[1]
    b_file = sys.argv[2]

    matrix = []
    with open(matrix_file) as m:
        for line in m:
            matrix.append(list(map(int, line.split())))

    Q, R = QR_decompostion(matrix)
    print(Q, R)
