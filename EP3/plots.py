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
    sigmas = [0.05, 0.1, 0.25, 0.5, 1, 5]
    for sigma in sigmas:
        Q0 = 2000
        L = 1
        Q1 = lambda x: Q0 * np.e ** (-((x - L / 2) ** 2) / sigma**2)
        Q1_lista = [Q1(x) for x in np.linspace(0, 1, 1000)]
        plt.plot(
            np.linspace(0, 1, 1000),
            Q1_lista,
            label="$\sigma$ = " + str(sigma),
        )
    plt.grid()
    plt.ylabel("Q(x)")
    plt.xlabel("x")
    plt.legend(loc="lower right")
    plt.show()

    thetas = [0.05, 0.1, 0.25, 0.5, 1, 5]
    for theta in thetas:
        Q_0 = 2000
        L = 1
        Q1 = lambda x: -(
            Q_0
            * (
                np.e ** (-(x**2) / theta**2)
                + np.e ** (-((x - L) ** 2) / theta**2)
            )
        )
        Q1_lista1 = [Q1(x) for x in np.linspace(0, 1, 1000)]
        plt.plot(
            np.linspace(0, 1, 1000),
            Q1_lista1,
            label="$\u03B8=$ " + str(theta),
        )
    plt.grid()
    plt.legend(loc="lower right")
    plt.ylabel("Q(x)")
    plt.xlabel("x")
    plt.show()

    sigmas = [0.1, 5]
    thetas = [0.1, 5]
    for theta in thetas:
        for sigma in sigmas:
            Q0 = 2000
            Q_0 = 1000
            L = 1
            Q1 = lambda x: Q0 * np.e ** (-((x - L / 2) ** 2) / sigma**2) - (
                Q_0
                * (
                    np.e ** (-(x**2) / theta**2)
                    + np.e ** (-((x - L) ** 2) / theta**2)
                )
            )
            Q1_lista2 = [Q1(x) for x in np.linspace(0, 1, 1000)]
            plt.plot(
                np.linspace(0, 1, 1000),
                Q1_lista2,
                label="$\u03B8=$ " + str(theta) + " $\sigma$ = " + str(sigma),
            )
    plt.grid()
    plt.ylabel("Q(x)")
    plt.xlabel("x")
    plt.legend(loc="lower right")
    plt.show()


if __name__ == "__main__":
    main()
