import numpy as np

from funcs import (
    LU,
    decompoeMatrizTridiag,
    sistemaTridiagLU,
    criaMatrizTridiagCircular,
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
        A = input("Matriz a ser decomposta: ")

        U, L = LU(A)

        print("Matriz A: ")
        print("Matriz L: ")
        print("Matriz U: ")

    else:
        raise Exception("Modo não existente")


if __name__ == "__main__":
    main()
