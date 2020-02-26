import sys
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


def transpose_matrix(A):
    n = len(A)
    m = len(A[0])
    result = [[0 for _ in range(n)] for _ in range(m)]
    for i in range(n):
        for j in range(m):
            result[j][i] = A[i][j]
    return result


def find_max(A):
    n = len(A)
    i_r, j_r = 0, 0
    maxx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j]) > maxx:
                maxx = abs(A[i][j])
                i_r = i
                j_r = j
    return i_r, j_r


def spin_method(A):
    n = len(A)
    A_new = A.copy()
    eigenvectors = [[0 if i != j else 1 for j in range(n)] for i in range(n)]
    while True:
        i_maxx, j_maxx = find_max(A_new)
        if A_new[i_maxx][i_maxx] == A_new[j_maxx][j_maxx]:
            phi = pi / 4
        else:
            phi = 0.5 * atan((2 * A_new[i_maxx][j_maxx]) / (A_new[i_maxx][i_maxx] - A_new[j_maxx][j_maxx]))

        U = [[0 if i != j else 1 for j in range(n)] for i in range(n)]
        U[i_maxx][i_maxx] = cos(phi)
        U[j_maxx][j_maxx] = cos(phi)
        U[i_maxx][j_maxx] = -sin(phi)
        U[j_maxx][i_maxx] = sin(phi)

        eigenvectors = matrix_product(eigenvectors, U)

        UT = transpose_matrix(U)
        A_new = matrix_product(matrix_product(UT, A_new), U)
        epsilon_k = 0
        for i in range(n):
            for j in range(i):
                epsilon_k += A_new[i][j] ** 2
        epsilon_k = sqrt(epsilon_k)
        if epsilon_k < epsilon:
            break
    eigenvalues = [round(A_new[i][i], 3) for i in range(n)]
    eigenvectors = [[round(eigenvectors[i][j], 3) for j in range(n)] for i in range(n)]
    return eigenvalues, eigenvectors


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

    print(spin_method(matrix))
