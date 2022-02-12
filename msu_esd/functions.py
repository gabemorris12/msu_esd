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
