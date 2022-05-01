import numpy as np

from funcs import (
    LU,
    sistemaTridiagLU,
    sistemaTridiagCiclico,
)


def main():
    print("Escolha um modo:")
    print("Modo 1 - Faz decomposição LU de uma matriz nxn")
    print("Modo 2 - Resolve um sistema tridiagonal")
    print("Modo 3 - Resolve sistema tridiagonal ciclico com matriz pré definida")
    print("Modo 4 - Resolve sistema tridiagonal ciclico com matriz escolhida")
    modo = int(input("Modo: "))

    if modo == 1:
        A = input("Matriz nxn a ser decomposta: ")

        U, L = LU(A)

        print("Matriz A: ", A)
        print("Matriz L: ", L)
        print("Matriz U: ", U)

    elif modo == 2:
        print("Sistema A*x = y")
        A = input("Matriz tridiagonal: ")
        y = input("Vetor y: ")

        x = sistemaTridiagLU(A, y)

        print("Solução: ", x)

    elif modo == 3:
        n = int(input("Dimensão da matriz tridiagonal: "))

        # Inicializa listas
        a = np.zeros(n, dtype=np.float64)
        b = np.zeros(n, dtype=np.float64)
        c = np.zeros(n, dtype=np.float64)
        d = np.zeros(n, dtype=np.float64)

        # Calcula a, b, c e d
        for i in range(n):
            a[i] = (
                (2 * (i + 1) - 1) / (4 * (i + 1))
                if (i + 1) != n
                else (2 * (i + 1) - 1) / (2 * (i + 1))
            )
            c[i] = 1 - a[i]
            b[i] = 2
            d[i] = np.cos((2 * np.pi * (i + 1) ** 2) / (n**2))

        x = sistemaTridiagCiclico(a, b, c, d)

        print("Solução: ", x)

    elif modo == 4:
        print("Insira os vetores que definem a matriz")
        a = input("vetor a: ")
        b = input("vetor b: ")
        c = input("vetor c: ")

        print("Insira o vetor d para o sistema A*x=d")
        d = input("vetor d: ")

        x = sistemaTridiagCiclico(a, b, c, d)

        print("Solução: ", x)
    else:
        raise Exception("Modo não existente")


if __name__ == "__main__":
    main()
