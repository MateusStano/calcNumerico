import numpy as np

# Resolve sistema tridiagonais
def sistemaTridiagLU(a, b, c, d, n):
    """
    a: vetor de coeficientes da diagonal secundaria inferior
    b: vetor de coeficientes da diagonal principal
    c: vetor de coeficientes da diagonal secundaria superior
    d: vetor resultado
    n: dimensão do sistema nxn

    """
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


# Calculo da integral simples ou dupla
def integral_simples_ou_dupla(
    f, a, b, pontos_e_pesos_x, c=lambda x: 0, d=lambda x: 1, pontos_e_pesos_y=None
):
    """
    f: função de x e y a ser integrada
    a: limite inferior do intervalo de integração em x
    b: limite superior do intervalo de integração em x
    pontos_e_pesos_x: pontos e pesos para o método de Gauss em x (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)
    c: limite inferior do intervalo de integração em y (função de x)
    d: limite superior do intervalo de integração em y (função de x)
    pontos_e_pesos_y: pontos e pesos para o método de Gauss em y (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)

    """

    pontos_e_pesos_y = (
        pontos_e_pesos_x if pontos_e_pesos_y is None else pontos_e_pesos_y
    )

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


def calcula_produto_interno_phis_diagonal_principal(
    k, q, past_xi, current_xi, next_xi, h, pontos_e_pesos_phi
):
    """
    k: função k(x)
    q: função q(x)
    past_xi: xi[i-1], xi anterior
    current_xi: xi[i], xi atual
    next_xi: xi[i+1], xi posterior
    h: tamanho do intervalo
    pontos_e_pesos_phi: pontos e pesos para o método de Gauss em phi (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)

    """

    # Calcula o produto interno de phii e phii
    phii_phii_1 = lambda x, y: k(x, y) + (x - past_xi) ** 2 * q(x, y)
    phii_phii_2 = lambda x, y: k(x, y) + (next_xi - x) ** 2 * q(x, y)
    return (1 / h) ** 2 * (
        integral_simples_ou_dupla(
            phii_phii_1, past_xi, current_xi, pontos_e_pesos_phi
        )  # intervalo de integração [past_xi, current_xi]
        + integral_simples_ou_dupla(
            phii_phii_2, current_xi, next_xi, pontos_e_pesos_phi
        )  # intervalo de integração [current_xi, next_xi]
    )


def calcula_produto_interno_phis_diagonal_secundaria_inferior(
    k, q, past_xi, current_xi, h, pontos_e_pesos_phi
):
    """
    k: função k(x)
    q: função q(x)
    past_xi: xi[i-1], xi anterior
    current_xi: xi[i], xi atual
    h: tamanho do intervalo
    pontos_e_pesos_phi: pontos e pesos para o método de Gauss em phi (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)
    """

    # Calcula o produto interno de phii e phii+1
    phii_phij = lambda x, y: k(x, y) + (current_xi - x) * (x - past_xi) * q(x, y)
    return -((1 / h) ** 2) * integral_simples_ou_dupla(
        phii_phij, past_xi, current_xi, pontos_e_pesos_phi
    )  # intervalo de integração [past_xi, current_xi]


def calcula_produto_interno_phis_diagonal_secundaria_superior(
    k, q, current_xi, next_xi, h, pontos_e_pesos_phi
):
    """
    k: função k(x)
    q: função q(x)
    current_xi: xi[i], xi atual
    next_xi: xi[i+1], xi posterior
    h: tamanho do intervalo
    pontos_e_pesos_phi: pontos e pesos para o método de Gauss em phi (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)
    """

    # Calcula o produto interno de phii e phii-1
    phij_phii = lambda x, y: k(x, y) + (next_xi - x) * (x - current_xi) * q(x, y)
    return -((1 / h) ** 2) * integral_simples_ou_dupla(
        phij_phii, current_xi, next_xi, pontos_e_pesos_phi
    )  # intervalo de integração [current_xi, next_xi]


def calcula_produto_interno_f_phi(
    f, past_xi, current_xi, next_xi, h, pontos_e_pesos_phi
):
    """
    f: função f(x)
    past_xi: xi[i-1], xi anterior
    current_xi: xi[i], xi atual
    next_xi: xi[i+1], xi posterior
    h: tamanho do intervalo
    pontos_e_pesos_phi: pontos e pesos para o método de Gauss em phi (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)
    """

    # Calcula o produto interno de f(x) e phii
    f_phii1 = lambda x, y: (x - past_xi) * f(x, y)
    f_phii2 = lambda x, y: (next_xi - x) * f(x, y)
    return (
        1
        / h
        * (
            integral_simples_ou_dupla(f_phii1, past_xi, current_xi, pontos_e_pesos_phi)
            + integral_simples_ou_dupla(
                f_phii2, current_xi, next_xi, pontos_e_pesos_phi
            )
        )
    )


def sol_sistema_linear_tridiagonal(f, k, q, n, pontos_e_pesos_phi, L=1):
    """
    f: função f(x)
    k: função k(x)
    q: função q(x)
    n: número de pontos
    pontos_e_pesos_phi: pontos e pesos para o método de Gauss em phi (deve ser da forma [(ponto1, peso1), (ponto2, peso2), ...].)
    L: opcional, define interval [0,L]
    """

    h = L / (n + 1)
    xi = [i * h for i in range(n + 2)]

    # Vetores do sistema linear (tridigonal)
    b = []  # diagonal principal
    a = []  # diagonal secundária inferior
    c = []  # diagonal secundária superior
    d = []  # vetor de termos independentes
    for i in range(n):
        past_xi = xi[i]
        current_xi = xi[i + 1]
        next_xi = xi[i + 2]

        b.append(
            calcula_produto_interno_phis_diagonal_principal(
                k, q, past_xi, current_xi, next_xi, h, pontos_e_pesos_phi
            )
        )
        a.append(
            calcula_produto_interno_phis_diagonal_secundaria_inferior(
                k, q, past_xi, current_xi, h, pontos_e_pesos_phi
            )
        )
        c.append(
            calcula_produto_interno_phis_diagonal_secundaria_superior(
                k, q, current_xi, next_xi, h, pontos_e_pesos_phi
            )
        )
        d.append(
            calcula_produto_interno_f_phi(
                f, past_xi, current_xi, next_xi, h, pontos_e_pesos_phi
            )
        )

    # Resolve o sistema linear
    alphas = np.array(sistemaTridiagLU(a, b, c, d, n))

    return alphas


def phi(x, xi, i, h):
    """
    x: valor de x
    xi: vetor de pontos
    i: índice do ponto de xi
    h: tamanho do intervalo
    """

    if xi[i - 1] < x <= xi[i]:
        return (x - xi[i - 1]) / h
    elif xi[i] < x <= xi[i + 1]:
        return (xi[i + 1] - x) / h
    else:
        return 0


def u_barra(x, alphas, xi, h):
    """
    x: valor de x
    alphas: vetor de alphas
    xi: vetor de pontos
    h: tamanho do intervalo

    """

    u_barra = 0
    for i in range(len(alphas)):
        u_barra += alphas[i] * phi(x, xi, i + 1, h)
    return u_barra


def u_barra_nao_homogenea(x, a, b, alphas, xi, h, L=1):
    """
    x: valor de x
    alphas: vetor de alphas
    xi: vetor de pontos
    h: tamanho do intervalo

    """
    return u_barra(x, alphas, xi, h) + (a + (b - a) * x) / L


def maior_erro(n, u_exato, alphas, L=1):
    """
    n: número de pontos
    u_exato: função u(x)
    alphas: vetor de alphas
    L: opcional, define interval [0,L]
    """

    h = L / (n + 1)
    xi = [i * h for i in range(n + 2)]
    erro_max = 0
    u_barra_erro_max = 0
    u_exato_erro_max = 0
    for x in xi:
        if abs(u_exato(x) - u_barra(x, alphas, xi, h)) > erro_max:
            erro_max = abs(u_exato(x) - u_barra(x, alphas, xi, h))
            u_barra_erro_max = u_barra(x, alphas, xi, h)
            u_exato_erro_max = u_exato(x)
    return erro_max, u_barra_erro_max, u_exato_erro_max
