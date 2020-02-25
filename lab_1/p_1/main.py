import sys


class MatrixError(Exception):
    pass


class InverseMatrixError(MatrixError):
    pass


class DegenerateMatrixError(MatrixError):
    pass


def get_cofactor(A, row, column, m):
    n = len(A)
    i1 = 0
    j1 = 0
    temp = [[0 for _ in range(m)] for _ in range(m)]
    for i in range(n):
        for j in range(n):
            if i != row and j != column:
                temp[i1][j1] = A[i][j]
                j1 += 1
                if j1 == m:
                    j1 = 0
                    i1 += 1
    return temp


def adjoint(A):
    n = len(A)
    adjoint_A = [[0 for _ in range(n)] for _ in range(n)]
    if n == 1:
        adjoint_A[0][0] = 1
        return adjoint_A

    sign = 1
    for i in range(n):
        for j in range(n):
            temp = get_cofactor(A, i, j, n - 1)

            sign = 1 if (i + j) % 2 == 0 else -1

            adjoint_A[j][i] = sign * determinant_of_matrix(temp)

    return adjoint_A


def inverse(A):
    n = len(A)
    det = determinant_of_matrix(A)
    if det == 0:
        raise InverseMatrixError("The inverse matrix does not exist")

    adjoint_A = adjoint(A)

    inverse_A = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            inverse_A[i][j] = adjoint_A[i][j] / det

    return inverse_A


def determinant_of_matrix(A):
    n = len(A)

    if n == 1:
        return A[0][0]

    sign = 1
    D = 0
    for i in range(n):
        temp = get_cofactor(A, 0, i, n - 1)
        D += sign * A[0][i] * determinant_of_matrix(temp)
        sign = -sign

    return D


def LUP_decomposition(A):
    n = len(A)
    pi = list(range(0, len(A)))
    for k in range(n):
        p = -1000
        for i in range(k, n):
            if A[i][k] > p:
                p = A[i][k]
                k_ = i
        if p == -1000:
            raise DegenerateMatrixError("Matrix is degenerate")
        pi[k], pi[k_] = pi[k_], pi[k]
        for i in range(n):
            A[k][i], A[k_][i] = A[k_][i], A[k][i]
        for i in range(k + 1, n):
            A[i][k] = A[i][k] / A[k][k]
            for j in range(k + 1, n):
                A[i][j] = A[i][j] - A[i][k] * A[k][j]

    L = list()
    U = list()
    for i in range(len(A)):
        L_ = []
        U_ = []
        for j in range(len(A)):
            if i > j:
                L_.append(round(A[i][j], 1))
                U_.append(0)
            else:
                if i == j:
                    L_.append(1)
                else:
                    L_.append(0)
                U_.append(round(A[i][j], 1))
        L.append(L_)
        U.append(U_)
    return L, U, pi


def LUP_solve(L, U, pi, b):
    n = len(L)
    x, y = [0 for _ in range(n)], [0 for _ in range(n)]
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += L[i][j] * y[j]
        y[i] = b[pi[i]] - sum

    for i in range(n - 1, -1, -1):
        sum = 0
        for j in range(i + 1, n):
            sum += U[i][j] * x[j]
        x[i] = round((y[i] - sum) / U[i][i], 1)
    return x


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("use {} <matrix_file> <b_file>")
        exit(0)

    matrix_file = sys.argv[1]
    b_file = sys.argv[2]

    A = []
    with open(matrix_file) as m:
        for line in m:
            A.append(list(map(int, line.split())))

    b = []
    with open(b_file) as m:
        b = list(map(int, m.read().split()))

    print("determinant A = {}".format(determinant_of_matrix(A)))
    print("inverse A:")
    inverse_A = inverse(A)
    for i in inverse_A:
        print(i)

    L, U, pi = LUP_decomposition(A)
    x = LUP_solve(L, U, pi, b)
    print("L:")
    for i in L:
        print(i)
    print("U:")
    for i in U:
        print(i)
    print("P:")
    for i in pi:
        for j in range(len(pi)):
            if i == j:
                print(1, end=" ")
            else:
                print(0, end=" ")
        print()
    print("x:")
    for i in x:
        print(i)
