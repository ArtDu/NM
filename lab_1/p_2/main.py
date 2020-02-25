import sys


class MatrixError(Exception):
    pass


def forward(matrix, D):
    n = len(matrix)
    A = [0 for _ in range(n)]
    B = [0 for _ in range(n)]
    b = matrix[0][0]
    c = matrix[0][1]
    d = D[0]
    A[0] = -c / b
    B[0] = d / b
    for i in range(1, n - 1):
        a = matrix[i][i - 1]
        b = matrix[i][i]
        c = matrix[i][i + 1]
        d = D[i]
        A[i] = -c / (b + a * A[i - 1])
        B[i] = (d - a * B[i - 1]) / (b + a * A[i - 1])
    A[n - 1] = 0
    a = matrix[n - 1][n - 2]
    b = matrix[n - 1][n - 1]
    d = D[n - 1]
    B[n - 1] = (d - a * B[n - 2]) / (b + a * A[n - 2])
    return A, B


def back(A, B):
    n = len(A)
    x = [0 for _ in range(n)]
    x[n - 1] = B[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = A[i] * x[i + 1] + B[i]
    return x


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

    b = []
    with open(b_file) as m:
        b = list(map(int, m.read().split()))

    A, B = forward(matrix, b)
    x = list(map(int, back(A, B)))
    print(x)
