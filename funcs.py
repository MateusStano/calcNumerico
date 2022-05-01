import numpy as np

# Decompõe uma matriz quadrada de dimensão n na for A = LU
def LU(A):
    # Dimensão n da matriz A
    n = A.shape[0]

    # Inicializa matrizes L e U
    L = np.identity(n, dtype=np.float64)
    U = np.zeros((n, n), dtype=np.float64)

    # Calcula U e L
    for i in range(n):
        U[i, i:] = A[i, i:] - L[i, :i] @ U[:i, i:]
        L[(i + 1) :, i] = (1 / U[i, i]) * (A[(i + 1) :, i] - L[(i + 1) :, :] @ U[:, i])

    return U, L


# Decompõe uma matriz tridiagonal quadrada de dimensão n na for A = LU
def decompoeMatrizTridiag(A):
    # Dimensão n da matriz A
    n = A.shape[0]

    # Define os vetores das diagonais

    # Diagonal secundaria abaixo da principal
    a = [A[i + 1][i] for i in range(n - 1)]
    a = [0] + a  # Adiciona 0 a primeiro item do vetor

    # Diagonal principal
    b = [A[i][i] for i in range(n)]

    # Diagonal secundaria acima da principal
    c = [A[i][i + 1] for i in range(n - 1)]
    c = c + [0]  # Adiciona 0 ao fim do vetor

    return a, b, c


# Resolve sistema tridiagonais
def sistemaTridiagLU(A, d):
    # Dimensão de A
    n = A.shape[0]

    # Define diagonais a, b e c
    a, b, c = decompoeMatrizTridiag(A)

    # Calcula vetores u e l
    u = [b[0]]
    l = []
    for i in range(1, n):
        l.append(a[i] / u[i - 1])
        u.append(b[i] - l[i - 1] * c[i - 1])

    # Calcula solução de L*y = d
    y = [d[0]]
    for i in range(1, n):
        y.append(d[i] - l[i - 1] * y[i - 1])

    # Calcula solução de U*x = y
    x = [0] * n
    x[n - 1] = y[n - 1] / u[n - 1]
    for i in reversed(range(0, n - 1)):
        x[i] = (y[i] - c[i] * x[i + 1]) / u[i]

    return x


# Cria matriz tridiagonal circular a partir dos vetores que definem suas diagonais
def criaMatrizTridiagCircular(a, b, c):
    # Dimensão n da matriz
    n = len(b)

    # Inicializa matriz de dimensão n
    A = np.zeros((n, n), dtype=np.float64)

    # Constrõe a matriz a partir de a,b e c
    for i in range(n):
        if i == 0:
            A[i][n - 1] = a[i]
        else:
            A[i][i - 1] = a[i]
        A[i][i] = b[i]
        if i == n - 1:
            A[i][0] = c[i]
        else:
            A[i][i + 1] = c[i]

    return A


# Resolve sistemas tridiagonais ciclicos
def sistemaTridiagCiclico(a, b, c, d):
    n = len(b)

    # Constroi matriz A a partir de a,b e c
    A = criaMatrizTridiagCircular(a, b, c)

    # Constroi submatriz principal T
    T = np.delete(A, n - 1, 1)
    T = np.delete(T, n - 1, 0)

    # Constroi v
    v = [0] * n
    v[0] = a[0]
    v[-1] = c[n - 2]

    # Resolve sistema Ty=d
    y = sistemaTridiagLU(T, d)

    # Resolve sistema Tz=v
    z = sistemaTridiagLU(T, v)

    # Solução do sistema
    x = [0] * n
    for i in range(n):
        x[i] = (d[i] - c[i] * y[0] - a[i] * y[i - 1]) / (
            b[i] - c[i] * z[0] - a[i] * z[i - 1]
        )

    return x, A
