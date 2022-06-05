"""
EP1 - Decomposição LU para Matrizes Tridiagonais

Nome: Laura do Prado Gonçalves Pinto
NUSP: 11819960

Nome: Mateus Stano Junqueira
NUSP: 11804845
"""

import numpy as np

from functions import (
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
        n = int(input("Dimensão da Matriz: "))
        A = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            for j in range(n):
                a_ij = float(input("item a" + str(i) + str(j) + ": "))
                A[i][j] = a_ij

        print("Matriz de entrada: ")
        print(A)

        U, L = LU(A)

        print("Matriz L: ")
        print(L)
        print("Matriz U: ")
        print(U)

    elif modo == 2:
        n = int(input("Dimensão da Matriz do sistema: "))
        A = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            for j in range(n):
                a_ij = float(input("Item a" + str(i) + str(j) + ": "))
                A[i][j] = a_ij

        print("Matriz de entrada: ")
        print(A)

        print("Insira itens do vetor y")
        y = [0] * n
        for i in range(n):
            yi = float(input("y" + str(i) + ": "))
            y[i] = yi

        print("Vetor y: ")
        print(y)

        x = sistemaTridiagLU(A, y)

        print("Solução x: ")
        print(x)

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

        x, A = sistemaTridiagCiclico(a, b, c, d)

        print("\nVetor a: ")
        print(a)
        print("\nVetor b: ")
        print(b)
        print("\nVetor c: ")
        print(c)
        print("\nVetor d: ")
        print(d)
        print("\nSolução x: ")
        for i in x:
            print(i)

    elif modo == 4:
        n = int(input("Dimensão da matriz: "))

        print("Insira itens do vetor a")
        a = [0] * n
        for i in range(n):
            ai = float(input("a" + str(i) + ": "))
            a[i] = ai

        print("\nVetor a: ")
        print(a)

        print("\nInsira itens do vetor b")
        b = [0] * n
        for i in range(n):
            bi = float(input("b" + str(i) + ": "))
            b[i] = bi

        print("\nVetor b: ")
        print(b)

        print("\nInsira itens do vetor c")
        c = [0] * n
        for i in range(n):
            ci = float(input("c" + str(i) + ": "))
            c[i] = ci

        print("\nVetor c: ")
        print(c)

        print("\nInsira itens do vetor d")
        d = [0] * n
        for i in range(n):
            di = float(input("d" + str(i) + ": "))
            d[i] = di

        print("\nVetor d: ")
        print(d)

        x, A = sistemaTridiagCiclico(a, b, c, d)

        print("\nSolução x: ")
        for i in x:
            print(i)

    else:
        raise Exception("Modo não existente")


if __name__ == "__main__":
    main()
