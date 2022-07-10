"""
EP3 - Modelagem de um Sistema de Resfriamento de Chips

Nome: Laura do Prado Gonçalves Pinto
NUSP: 11819960

Nome: Mateus Stano Junqueira
NUSP: 11804845
"""

from functions import *
import numpy as np
import matplotlib.pyplot as plt

# pontos e pesos pré estabelecidos
# na forma de n = [(x_j,w_j),...]
n6 = [
    (0.2386191860831969086305017, 0.4679139345726910473898703),
    (0.6612093864662645136613996, 0.3607615730481386075698335),
    (0.9324695142031520278123016, 0.1713244923791703450402961),
    (-0.2386191860831969086305017, 0.4679139345726910473898703),
    (-0.6612093864662645136613996, 0.3607615730481386075698335),
    (-0.9324695142031520278123016, 0.1713244923791703450402961),
]
n8 = [
    (0.1834346424956498049394761, 0.3626837833783619829651504),
    (0.5255324099163289858177390, 0.3137066458778872873379622),
    (0.7966664774136267395915539, 0.2223810344533744705443560),
    (0.9602898564975362316835609, 0.1012285362903762591525314),
    (-0.1834346424956498049394761, 0.3626837833783619829651504),
    (-0.5255324099163289858177390, 0.3137066458778872873379622),
    (-0.7966664774136267395915539, 0.2223810344533744705443560),
    (-0.9602898564975362316835609, 0.1012285362903762591525314),
]
n10 = [
    (0.1488743389816312108848260, 0.2955242247147528701738930),
    (0.4333953941292471907992659, 0.2692667193099963550912269),
    (0.6794095682990244062343274, 0.2190863625159820439955349),
    (0.8650633666889845107320967, 0.1494513491505805931457763),
    (0.9739065285171717200779640, 0.0666713443086881375935688),
    (-0.1488743389816312108848260, 0.2955242247147528701738930),
    (-0.4333953941292471907992659, 0.2692667193099963550912269),
    (-0.6794095682990244062343274, 0.2190863625159820439955349),
    (-0.8650633666889845107320967, 0.1494513491505805931457763),
    (-0.9739065285171717200779640, 0.0666713443086881375935688),
]


def main():
    print("Escolha um modo:")
    print("Modo 1 - Validação no intervalo [0,1] com k(x)=1, q(x)=0, f(x)=12x(1-x)-2")
    print("Modo 2 - Validação no intervalo [0,1] com k(x)=e**x, q(x)=0, f(x)=e**x+1")
    print("Modo 3 - Equilíbrio com forçantes de calor")
    print("Modo 4 - Resolve sistema tridiagonal ciclico com matriz escolhida")
    modo = int(input("Modo: "))

    if modo == 1:
        print(
            "\nModo 1 - Validação no intervalo [0,1] com k(x)=1, q(x)=0, f(x)=12x(1-x)-2"
        )
        print("\nSolução exata u(x) = (x**2 *(1 - x)**2")
        print("\nParametros: \nk(x)=1, \nq(x)=0, \nf(x)=12x(1-x)-2")

        print("\n\nAvaliando o máximo erro para diferentes n's:")

        for n in [7, 15, 31, 63]:  # números de pontos

            print(f"\nn={n} -->")

            f = lambda x, y: 12 * x * (1 - x) - 2  # função f(x)
            k = lambda x, y: 1  # função k(x)
            q = lambda x, y: 0  # função q(x)
            u_exato = lambda x: x**2 * (1 - x) ** 2  # solução exata
            alphas = sol_sistema_linear_tridiagonal(
                f, k, q, n, n10
            )  # solução do sistema

            erro_max, u_barra_erro_max, u_exato_erro_max = maior_erro(
                n, u_exato, alphas
            )
            print(
                f"u_barra = {u_barra_erro_max:.22f} ; u_exato = {u_exato_erro_max:.22f}; erro = {erro_max:.22f}"
            )

        input("\nTodos os erros calculados. Pressione Enter para continuar...")
        for n in [2, 3, 7, 15, 31, 63]:

            print(f"\nn={n} -->")

            alphas = sol_sistema_linear_tridiagonal(
                f, k, q, n, n10
            )  # solução do sistema
            L = 1
            h = L / (n + 1)
            xi = [i * h for i in range(n + 2)]
            u_barra_lista = []
            for x in xi:
                u_barra_lista.append(u_barra(x, alphas, xi, h))
                print(
                    f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f} ; u_exato = {u_exato(x):.22f}; erro = {u_barra(x,alphas,xi,h)-u_exato(x):.22f}"
                )
            u_exato_lista = [u_exato(x) for x in np.linspace(0, 1)]

            plt.plot(xi, u_barra_lista, label="u_barra")
            plt.plot(np.linspace(0, 1), u_exato_lista, label="u_exato")
            plt.show()
    elif modo == 2:
        print(
            "\nModo 1 - Validação no intervalo [0,1] com k(x)=e**x, q(x)=0, f(x)=e**(x)+1"
        )
        print("\nSolução exata u(x) = (x-1)*(e**(-x) - 1)")
        print("\nParametros: \nk(x)=1, \nq(x)=0, \nf(x)=12x(1-x)-2")

        print("\n\nAvaliando o máximo erro para diferentes n's:")

        for n in [7, 15, 31, 63]:  # números de pontos

            print(f"\nn={n} -->")

            f = lambda x, y: np.exp(x) + 1  # função f(x)
            k = k = lambda x, y: np.exp(x)  # função k(x)
            q = lambda x, y: 0  # função q(x)
            u_exato = lambda x: (x - 1) * (np.exp(-x) - 1)  # solução exata
            alphas = sol_sistema_linear_tridiagonal(
                f, k, q, n, n10
            )  # solução do sistema

            erro_max, u_barra_erro_max, u_exato_erro_max = maior_erro(
                n, u_exato, alphas
            )
            print(
                f"u_barra = {u_barra_erro_max:.22f} ; u_exato = {u_exato_erro_max:.22f}; erro = {erro_max:.22f}"
            )

        input("\nTodos os erros calculados. Pressione Enter para continuar...")
        for n in [7, 15, 31, 63]:

            print(f"\nn={n} -->")

            alphas = sol_sistema_linear_tridiagonal(
                f, k, q, n, n10
            )  # solução do sistema
            L = 1
            h = L / (n + 1)
            xi = [i * h for i in range(n + 2)]
            for x in xi:
                print(
                    f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f} ; u_exato = {u_exato(x):.22f}; erro = {u_barra(x,alphas,xi,h)-u_exato(x):.22f}"
                )

    elif modo == 3:
        print("\nModo 3 - Equilíbrio com forçantes de calor")
        print("\nParametros: \nk(x)=3.6, \nq(x)=0, \nf(x)=Q(x)")
        print(
            "\nSendo Q(x) representa a soma do calor gerado pelo chip (Q+) e o calor retirado pelo resfriador (Q-)"
        )
        print("\n\nVariando parametros:")

        print(
            "\n  Q+ e Q- constantes e Q+ > Q- (portanto Q(x) constante e maior que zero)"
        )
        print(" Parametros: Q+ - Q- = 2")

        n = 15
        f = lambda x, y: 2  # função f(x)
        k = lambda x, y: 3.6  # função k(x)
        q = lambda x, y: 0  # função q(x)
        alphas = sol_sistema_linear_tridiagonal(f, k, q, n, n10)  # solução do sistema

        L = 1
        h = L / (n + 1)
        xi = [i * h for i in range(n + 2)]
        u_barra_lista = []
        for x in xi:
            u_barra_lista.append(u_barra(x, alphas, xi, h))
            print(f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f}")

        plt.plot(xi, u_barra_lista, label="$T(x)$")
        plt.legend()
        plt.grid()
        plt.show()

        input("\n\n  Q+ e Q- constantes e Q+ >>> Q-. Pressione Enter para continuar...")
        print(" Parametros: Q+ - Q- = 20")

        n = 15
        f = lambda x, y: 20  # função f(x)
        k = lambda x, y: 3.6  # função k(x)
        q = lambda x, y: 0  # função q(x)
        alphas = sol_sistema_linear_tridiagonal(f, k, q, n, n10)  # solução do sistema

        L = 1
        h = L / (n + 1)
        xi = [i * h for i in range(n + 2)]
        u_barra_lista = []
        for x in xi:
            u_barra_lista.append(u_barra(x, alphas, xi, h))
            print(f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f}")

        plt.plot(xi, u_barra_lista, label="$T(x)$")
        plt.legend()
        plt.grid()
        plt.show()

        input(
            "\n\n  Q+ e Q- constantes e Q+ < Q- (portanto Q(x) constante e menor que zero). Pressione Enter para continuar..."
        )
        print(" Parametros: Q+ - Q- = -2")

        n = 15
        f = lambda x, y: -2  # função f(x)
        k = lambda x, y: 3.6  # função k(x)
        q = lambda x, y: 0  # função q(x)
        alphas = sol_sistema_linear_tridiagonal(f, k, q, n, n10)  # solução do sistema

        L = 1
        h = L / (n + 1)
        xi = [i * h for i in range(n + 2)]
        u_barra_lista = []
        for x in xi:
            u_barra_lista.append(u_barra(x, alphas, xi, h))
            print(f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f}")

        plt.plot(xi, u_barra_lista, label="$T(x)$")
        plt.legend()
        plt.grid()
        plt.show()

        input("\n\n  Q+ e Q- constantes e Q+ <<< Q-. Pressione Enter para continuar...")
        print(" Parametros: Q+ - Q- = -20")

        n = 15
        f = lambda x, y: -20  # função f(x)
        k = lambda x, y: 3.6  # função k(x)
        q = lambda x, y: 0  # função q(x)
        alphas = sol_sistema_linear_tridiagonal(f, k, q, n, n10)  # solução do sistema

        L = 1
        h = L / (n + 1)
        xi = [i * h for i in range(n + 2)]
        u_barra_lista = []
        for x in xi:
            u_barra_lista.append(u_barra(x, alphas, xi, h))
            print(f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f}")

        plt.plot(xi, u_barra_lista, label="$T(x)$")
        plt.legend()
        plt.grid()
        plt.show()

        input(
            "\n\n  Q+ modelado por Q+(x) = Q+0 * e**-(x-L/2)**2/σ**2 e Q- constante. Pressione Enter para continuar..."
        )
        print(" Parametros: Q+0 = 2, σ = 0.5, Q-0 = 2")
        n = 15
        Q0 = 2
        σ = 0.5
        Q_0 = 2
        f = (
            lambda x, y: Q0 * np.e ** (-((x - L / 2) ** 2) / σ**2) - Q_0
        )  # função f(x)
        k = lambda x, y: 3.6  # função k(x)
        q = lambda x, y: 0  # função q(x)
        alphas = sol_sistema_linear_tridiagonal(f, k, q, n, n10)  # solução do sistema

        L = 1
        h = L / (n + 1)
        xi = [i * h for i in range(n + 2)]
        u_barra_lista = []
        for x in xi:
            u_barra_lista.append(u_barra(x, alphas, xi, h))
            print(f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f}")

        plt.plot(xi, u_barra_lista, label="$T(x)$")
        plt.legend()
        plt.grid()
        plt.show()

        input(
            "\n\n  Q+ constance e Q- modelado por Q-(x) = Q-0 * (e**(-x**2/θ**2) + e**(-(x-L)**2/θ/**2)). Pressione Enter para continuar..."
        )

        n = 15
        f = lambda x, y: -2  # função f(x)
        k = lambda x, y: 3.6  # função k(x)
        q = lambda x, y: 0  # função q(x)
        alphas = sol_sistema_linear_tridiagonal(f, k, q, n, n10)  # solução do sistema

        L = 1
        h = L / (n + 1)
        xi = [i * h for i in range(n + 2)]
        u_barra_lista = []
        for x in xi:
            u_barra_lista.append(u_barra(x, alphas, xi, h))
            print(f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f}")

        plt.plot(xi, u_barra_lista, label="$T(x)$")
        plt.legend()
        plt.grid()
        plt.show()

        input(
            "\n\n  Q+ modelado por Q+(x) = Q+0 * e**-(x-L/2)**2/σ**2 e Q- modelado por Q-(x) = Q-0 * (e**(-x**2/θ**2) + e**(-(x-L)**2/θ/**2)). Pressione Enter para continuar..."
        )

        n = 15
        f = lambda x, y: -2  # função f(x)
        k = lambda x, y: 3.6  # função k(x)
        q = lambda x, y: 0  # função q(x)
        alphas = sol_sistema_linear_tridiagonal(f, k, q, n, n10)  # solução do sistema

        L = 1
        h = L / (n + 1)
        xi = [i * h for i in range(n + 2)]
        u_barra_lista = []
        for x in xi:
            u_barra_lista.append(u_barra(x, alphas, xi, h))
            print(f"x = {x:.4f} ; u_barra = {u_barra(x,alphas,xi,h):.22f}")

        plt.plot(xi, u_barra_lista, label="$T(x)$")
        plt.legend()
        plt.grid()
        plt.show()


if __name__ == "__main__":
    main()
