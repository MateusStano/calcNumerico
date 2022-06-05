"""
EP2 - Formulas de Integracao Numerica de Gauss

Nome: Laura do Prado Gonçalves Pinto
NUSP: 11819960

Nome: Mateus Stano Junqueira
NUSP: 11804845
"""

import numpy as np

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


def integral_dupla(f, a, b, pontos_e_pesos_x, c, d, pontos_e_pesos_y):
    """
    f: função de x e y a ser integrada
    a: limite inferior do intervalo de integração em x
    b: limite superior do intervalo de integração em x
    pontos_e_pesos_x: pontos e pesos para o método de Gauss em x (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)
    c: limite inferior do intervalo de integração em y (função de x)
    d: limite superior do intervalo de integração em y (função de x)
    pontos_e_pesos_y: pontos e pesos para o método de Gauss em y (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)

    """

    I = 0
    for x_i, w_i in pontos_e_pesos_x:
        # troca de variavel de x para o intervalo [a,b]
        x_i = (a + b) / 2 + (b - a) * x_i / 2
        F = 0
        for y_ij, w_ij in pontos_e_pesos_y:
            # troca de variavel de y para o intervalo [c(x_i),d(x_i)]
            y_ij = (c(x_i) + d(x_i)) / 2 + (d(x_i) - c(x_i)) * y_ij / 2
            F += w_ij * f(x_i, y_ij) * (d(x_i) - c(x_i)) / 2
        I += w_i * F * (b - a) / 2

    return I


def main():
    print("Escolha um modo:")
    print("Modo 0 - Calculo de integral dupla para uma função e intervalos quaisquer")
    print("Modo 1 - Resolve o Exemplo 1")
    print("Modo 2 - Resolve o Exemplo 2")
    print("Modo 3 - Resolve o Exemplo 3")
    print("Modo 4 - Resolve o Exemplo 4")
    modo = int(input("Modo: "))

    if modo == 0:
        """Modo 0 - Calculo de integral dupla para uma função e intervalos quaisquer"""

        print("Calculo de integral dupla da forma dydx")
        f = input("Função f(x,y) a ser integrada: ")
        f_str = ""
        for char in f:  # Converte a função para uma string utilizavel em eval()
            if char == "^":
                f_str += "**"
            else:
                f_str += char
        f_final = lambda x, y: eval(
            f
        )  # Cria uma funcao lambda que recebe x e y e retorna o valor da função f(x,y)

        a = float(
            input(
                "Limite inferior do intervalo de integração da integral exterior (dx): "
            )
        )

        b = float(
            input(
                "Limite superior do intervalo de integração da integral exterior (dx): "
            )
        )

        nx = int(input("Número de subintervalos (nós) para x (opções 6, 8 ou 10): "))
        if nx == 6:
            nx = n6
        elif nx == 8:
            nx = n8
        elif nx == 10:
            nx = n10
        else:
            raise ValueError("Número de subintervalos inválido")

        c = input(
            "Limite inferior do intervalo de integração da integral exterior (dy): "
        )
        c_str = ""
        for char in c:  # Converte a função para uma string utilizavel em eval()
            if char == "^":
                c_str += "**"
            else:
                c_str += char
        c_final = lambda x: eval(
            c_str
        )  # Cria uma funcao lambda que recebe x e retorna o valor da função c(x)

        d = input(
            "Limite superior do intervalo de integração da integral exterior (dy): "
        )
        d_str = ""
        for char in d:  # Converte a função para uma string utilizavel em eval()
            if char == "^":
                d_str += "**"
            else:
                d_str += char
        d_final = lambda x: eval(
            d_str
        )  # Cria uma funcao lambda que recebe x e retorna o valor da função d(x)

        ny = int(input("Número de subintervalos (nós) para y (opções 6, 8 ou 10): "))
        if ny == 6:
            ny = n6
        elif ny == 8:
            ny = n8
        elif ny == 10:
            ny = n10
        else:
            raise ValueError("Número de subintervalos inválido")

        print("\nIntegral da forma dydx:")
        print(
            "I = ",
            integral_dupla(f_final, a, b, nx, c_final, d_final, ny),
            "\n",
        )

    elif modo == 1:
        """Modo 1 - Resolve o Exemplo 1"""

        # Cubo ------------------------------------------------------------------
        print("\n\nCálculo do volume do cubo cujas arestas tem comprimento 1:\n")
        print("Função a ser integrada: f(x,y) = 1\n")
        print("Intervalos de integração: \n a = 0, b = 1, c = 0, d = 1\n")
        print("Resultados: ")
        print(
            "n = 6 -->  I = ",
            integral_dupla(
                f=lambda x, y: 1,
                a=0,
                b=1,
                pontos_e_pesos_x=n6,
                c=lambda x: 0,
                d=lambda x: 1,
                pontos_e_pesos_y=n6,
            ),
        )
        print(
            "n = 8 -->  I = ",
            integral_dupla(
                f=lambda x, y: 1,
                a=0,
                b=1,
                pontos_e_pesos_x=n8,
                c=lambda x: 0,
                d=lambda x: 1,
                pontos_e_pesos_y=n8,
            ),
        )
        print(
            "n = 10 --> I = ",
            integral_dupla(
                f=lambda x, y: 1,
                a=0,
                b=1,
                pontos_e_pesos_x=n10,
                c=lambda x: 0,
                d=lambda x: 1,
                pontos_e_pesos_y=n10,
            ),
            "\n",
        )

        print("Resultado exato: 1\n")

        # Tetraedro ------------------------------------------------------------------
        print(
            "\n\nCálculo do volume do tetraedro com vertices (0, 0, 0), (1, 0, 0), (0, 1, 0) e (0, 0, 1):\n"
        )
        print("Função a ser integrada: f(x,y) = 1-x-y\n")
        print("Intervalos de integração: \n a = 0, b = 1, c = 0, d = 1-x\n")
        print("Resultados: ")
        print(
            "n = 6 -->  I = ",
            integral_dupla(
                f=lambda x, y: 1 - x - y,
                a=0,
                b=1,
                pontos_e_pesos_x=n6,
                c=lambda x: 0,
                d=lambda x: 1 - x,
                pontos_e_pesos_y=n6,
            ),
        )
        print(
            "n = 8 -->  I = ",
            integral_dupla(
                f=lambda x, y: 1 - x - y,
                a=0,
                b=1,
                pontos_e_pesos_x=n8,
                c=lambda x: 0,
                d=lambda x: 1 - x,
                pontos_e_pesos_y=n8,
            ),
        )
        print(
            "n = 10 --> I = ",
            integral_dupla(
                f=lambda x, y: 1 - x - y,
                a=0,
                b=1,
                pontos_e_pesos_x=n10,
                c=lambda x: 0,
                d=lambda x: 1 - x,
                pontos_e_pesos_y=n10,
            ),
            "\n",
        )

        print("Resultado exato: 1/6\n")

    elif modo == 2:
        """Modo 2 - Resolve o Exemplo 2"""

        print(
            "\n\nA area da regiao no primeiro quadrante limitada pelos eixos e pela curva y = 1-x^2:\n"
        )
        print("Função a ser integrada: f(x,y) = 1\n")
        print("Intervalos de integração: \n a = 0, b = 1, c = 0, d = 1-x^2\n")
        print("Resultados: ")
        print(
            "n = 6 -->  I = ",
            integral_dupla(
                f=lambda x, y: 1,
                a=0,
                b=1,
                pontos_e_pesos_x=n6,
                c=lambda x: 0,
                d=lambda x: 1 - x**2,
                pontos_e_pesos_y=n6,
            ),
        )
        print(
            "n = 8 -->  I = ",
            integral_dupla(
                f=lambda x, y: 1,
                a=0,
                b=1,
                pontos_e_pesos_x=n8,
                c=lambda x: 0,
                d=lambda x: 1 - x**2,
                pontos_e_pesos_y=n8,
            ),
        )
        print(
            "n = 10 --> I = ",
            integral_dupla(
                f=lambda x, y: 1,
                a=0,
                b=1,
                pontos_e_pesos_x=n10,
                c=lambda x: 0,
                d=lambda x: 1 - x**2,
                pontos_e_pesos_y=n10,
            ),
            "\n",
        )

        print("Resultado exato: 2/3\n")

    elif modo == 3:
        """Modo 3 - Resolve o Exemplo 3"""
        print("Area e volume da superfície descrita por:")
        print("z = e^(y/x),")
        print("0.1 <= x <= 0.5,")
        print("x^3 <= y <= x^2\n")

        print("Obs: x e y devem ser trocados para ordem correta de integração\n")

        # Superficie ------------------------------------------------------------------
        print("Cálculo da superfície descrita:\n")
        print("Função a ser integrada: f(x,y) = e^(y/x)\n")
        print("Intervalos de integração: \n a = 0.1, b = 0.5, c = x^3, d = x^2\n")
        print("Resultados: ")
        print(
            "n = 6 -->  I = ",
            integral_dupla(
                f=lambda x, y: np.e ** (y / x),
                a=0.1,
                b=0.5,
                pontos_e_pesos_x=n6,
                c=lambda x: x**3,
                d=lambda x: x**2,
                pontos_e_pesos_y=n6,
            ),
        )
        print(
            "n = 8 -->  I = ",
            integral_dupla(
                f=lambda x, y: np.e ** (y / x),
                a=0.1,
                b=0.5,
                pontos_e_pesos_x=n8,
                c=lambda x: x**3,
                d=lambda x: x**2,
                pontos_e_pesos_y=n8,
            ),
        )
        print(
            "n = 10 --> I = ",
            integral_dupla(
                f=lambda x, y: np.e ** (y / x),
                a=0.1,
                b=0.5,
                pontos_e_pesos_x=n10,
                c=lambda x: x**3,
                d=lambda x: x**2,
                pontos_e_pesos_y=n10,
            ),
            "\n",
        )

        print("Resultado exato: 1.6\n")

        # Volume ------------------------------------------------------------------
        print("Cálculo do volume da superfície descrita:\n")
        print(
            "Função a ser integrada: f(x,y) = ( (e^(x/y)/y)^2 + (-x*e^(x/y)/y^2)^2 + 1 )^(1/2)\n"
        )
        print("Intervalos de integração: \n a = 0.1, b = 0.5, c = x^3, d = x^2\n")
        print("Resultados: ")
        print(
            "n = 6 -->  I = ",
            integral_dupla(
                f=lambda x, y: (
                    (np.e ** (y / x) / x) ** 2
                    + (-y * np.e ** (y / x) / x**2) ** 2
                    + 1
                )
                ** (1 / 2),
                a=0.1,
                b=0.5,
                pontos_e_pesos_x=n6,
                c=lambda x: x**3,
                d=lambda x: x**2,
                pontos_e_pesos_y=n6,
            ),
        )
        print(
            "n = 8 -->  I = ",
            integral_dupla(
                f=lambda x, y: (
                    (np.e ** (y / x) / x) ** 2
                    + (-y * np.e ** (y / x) / x**2) ** 2
                    + 1
                )
                ** (1 / 2),
                a=0.1,
                b=0.5,
                pontos_e_pesos_x=n8,
                c=lambda x: x**3,
                d=lambda x: x**2,
                pontos_e_pesos_y=n8,
            ),
        )
        print(
            "n = 10 --> I = ",
            integral_dupla(
                f=lambda x, y: (
                    (np.e ** (y / x) / x) ** 2
                    + (-y * np.e ** (y / x) / x**2) ** 2
                    + 1
                )
                ** (1 / 2),
                a=0.1,
                b=0.5,
                pontos_e_pesos_x=n10,
                c=lambda x: x**3,
                d=lambda x: x**2,
                pontos_e_pesos_y=n10,
            ),
            "\n",
        )

    elif modo == 4:
        """Modo 4 - Resolve o Exemplo 4"""

        # Volume da calota esferica --------------------------------------------------
        print(
            "\n\nCalculo do volume da calota esferica de altura 1/4 da esfera de raio 1\n"
        )
        print("Função a ser integrada: f(x,y) = y\n")
        print("Intervalos de integração: \n a = 3/4, b = 1, c = 0, d = (1-x^2)^(1/2)\n")
        print("Resultados: ")
        print(
            "n = 6 -->  I = ",
            2
            * np.pi
            * integral_dupla(
                f=lambda x, y: y,
                a=3 / 4,
                b=1,
                pontos_e_pesos_x=n6,
                c=lambda x: 0,
                d=lambda x: (1 - x**2) ** (1 / 2),
                pontos_e_pesos_y=n6,
            ),
        )
        print(
            "n = 8 -->  I = ",
            2
            * np.pi
            * integral_dupla(
                f=lambda x, y: y,
                a=3 / 4,
                b=1,
                pontos_e_pesos_x=n8,
                c=lambda x: 0,
                d=lambda x: (1 - x**2) ** (1 / 2),
                pontos_e_pesos_y=n8,
            ),
        )
        print(
            "n = 10 --> I = ",
            2
            * np.pi
            * integral_dupla(
                f=lambda x, y: y,
                a=3 / 4,
                b=1,
                pontos_e_pesos_x=n10,
                c=lambda x: 0,
                d=lambda x: (1 - x**2) ** (1 / 2),
                pontos_e_pesos_y=n10,
            ),
            "\n",
        )

        # Volume do solido de revolucao --------------------------------------------------
        print(
            "\n\nCalculo do volume do solido de revolucao obtido da rotacao da regiao x=0, x=e^(-y^2), y=-1, y=1 em torno do eixo y\n"
        )
        print("Função a ser integrada: f(x,y) = y\n")
        print("Intervalos de integração: \n a = -1, b = 1, c = 0, d = e^(-y^2)\n")
        print("Resultados: ")
        print(
            "n = 6 -->  I = ",
            2
            * np.pi
            * integral_dupla(
                f=lambda x, y: y,
                a=-1,
                b=1,
                pontos_e_pesos_x=n6,
                c=lambda x: 0,
                d=lambda x: np.exp(-(x**2)),
                pontos_e_pesos_y=n6,
            ),
        )
        print(
            "n = 8 -->  I = ",
            2
            * np.pi
            * integral_dupla(
                f=lambda x, y: y,
                a=-1,
                b=1,
                pontos_e_pesos_x=n8,
                c=lambda x: 0,
                d=lambda x: np.exp(-(x**2)),
                pontos_e_pesos_y=n8,
            ),
        )
        print(
            "n = 10 --> I = ",
            2
            * np.pi
            * integral_dupla(
                f=lambda x, y: y,
                a=-1,
                b=1,
                pontos_e_pesos_x=n10,
                c=lambda x: 0,
                d=lambda x: np.exp(-(x**2)),
                pontos_e_pesos_y=n10,
            ),
            "\n",
        )


if __name__ == "__main__":
    main()
