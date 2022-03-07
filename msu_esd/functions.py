import numpy as np


def f_T(D, epsilon):
    """
    Calculates the friction factor as Reynold's number approaches infinity. In other words, this is the friction factor
    for a fully rough flow. This utilizes the Haaland equations.

    :param D: The diameter of the pipe
    :param epsilon: The absolute roughness of the pipe
    :return: f_T
    """
    return 0.3086/(np.log10((epsilon/(3.7*D))**1.11))**2


def Re(Q, D, rho, mu):
    """
    Calculates the Reynolds number in terms of the flow rate. Everything must be in the base units of the system.

    :param Q: The flow rate
    :param D: The diameter
    :param rho: The density of the fluid
    :param mu: The dynamic viscosity of the fluid
    :return: Re
    """
    return (4*rho*np.abs(Q))/(np.pi*D*mu)


def f(Q, D, epsilon, rho, mu):
    """
    Calculates the friction factor as a function of the flow rate and other typically known parameters. Only the flow
    rate can be an array, the other parameters must be constants. The friction factor is calculated using the Haaland
    equations.

    :param Q: The flow rate
    :param D: The diameter
    :param epsilon: The absolute roughness of the pipe
    :param rho: The density of the fluid
    :param mu: The dynamic viscosity of the fluid
    :return: f
    """
    Q = np.float64(Q)  # This calculation raises some issues if Q is a large python integer. This fixes it.

    def turbulent(D_, epsilon_, rho_, mu_):
        return lambda Q__: 0.3086/(np.log10(6.9/Re(Q__, D_, rho_, mu_) + (epsilon_/(3.7*D_))**1.11))**2

    def laminar(D_, rho_, mu_):
        return lambda Q__: 64/Re(Q__, D_, rho_, mu_)

    return np.piecewise(Q, [Re(Q, D, rho, mu) > 2300, Re(Q, D, rho, mu) <= 2300],
                        [turbulent(D, epsilon, rho, mu), laminar(D, rho, mu)])


def hardy_cross(pipes, Q, N, h=None, dh=None, tol=0.0001):
    """
    Returns the Hardy Cross solution for a pipe arrangement.

    :param pipes: A list of pipe objects
    :param Q: The corresponding guess flow rates to the pipe objects list. Must satisfy equilibrium at the nodes.
    :param N: The connection matrix. Should be a "number of pipes" by "number of loops" shape.
    :param h: A list of callable functions. Only necessary if adding a device that adds head loss or gain.
    :param dh: A list of callable functions of the derivative of h.
    :param tol: The error tolerance
    :return: Q, a list of flow rates
    """
    assert isinstance(N, np.ndarray), '"N" must be a numpy array.'

    if h is None and dh is None:
        h = [lambda Q_: 0 for _ in range(len(Q))]
        dh = [lambda Q_: 0 for _ in range(len(Q))]

    P, L = N.shape
    del_Q = np.full(L, 100)
    r = np.sqrt(np.sum(del_Q**2))

    get_h = lambda: np.array([pipe.h(Q[i]) + h[i](Q[i]) for i, pipe in enumerate(pipes)])
    get_dh = lambda: np.array([pipe.dh(Q[i]) + dh[i](Q[i]) for i, pipe in enumerate(pipes)])

    while r > tol:
        del_Q = -1*np.matmul(get_h(), N)/np.matmul(get_dh(), N**2)
        r = np.sqrt(np.sum(del_Q**2))
        Q = (Q.reshape((P, 1)) + np.matmul(N, del_Q.reshape((L, 1)))).reshape(P)

    return Q


def log_mean_temp_difference(T_h_in, T_h_out, T_c_in, T_c_out):
    """
    Finds the log mean temperature difference

    :param T_h_in: Temp of hot fluid as it goes in
    :param T_h_out: Temp of hot fluid as it goes out
    :param T_c_in: Temp of cold fluid as it goes in
    :param T_c_out: Temp of cold fluid as it goes out
    :return: The log mean temperature difference
    """
    T2 = T_h_out - T_c_in
    T1 = T_h_in - T_c_out
    return (T2 - T1)/np.log(T2/T1)
