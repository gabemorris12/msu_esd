from .functions import f, f_T, Re
import numpy as np


class Pipe:
    def __init__(self, D, L, epsilon, rho, mu, K=0, C=0, P_in=0, P_out=0, V_in=0, V_out=0, z_in=0, z_out=0,
                 g=32.174, eta=1, pump=False, turbine=False):
        self.D, self.L, self.epsilon, self.rho, self.mu = D, L, epsilon, rho, mu
        self.K, self.C = K, C
        self.P_in, self.P_out = P_in, P_out
        self.V_in, self.V_out = V_in, V_out
        self.z_in, self.z_out = z_in, z_out
        self.g = g
        self.eta = eta
        self.pump, self.turbine = pump, turbine
        self.W_s_turbine, self.W_s_pump, self.Q = None, None, None

        # todo: This needs to be tested by incorporating problems that have both pumps and turbines as well as problems
        # todo: where the pump and turbine heads are known
        # True means it exists and is unknown, false means that it's zero
        if pump is True and turbine is True:
            self.solve = lambda W_s_turbine, W_s_pump, Q: self._energy_balance(W_s_turbine, W_s_pump, Q)
            self.unknowns = 3
        elif pump is True and turbine is False:
            self.solve = lambda W_s_pump, Q: self._energy_balance(0, W_s_pump, Q)
            self.unknowns = 2
        elif pump is False and turbine is True:
            self.solve = lambda W_s_turbine, Q: self._energy_balance(W_s_turbine, 0, Q)
            self.unknowns = 2
        elif pump is False and turbine is False:
            self.solve = lambda Q: self._energy_balance(0, 0, Q)
            self.unknowns = 1
        else:
            # If we get here then that means that either the pump or the turbine is not a boolean. If it's not a
            # boolean, then I want them to be interpreted as the work per unit mass (W_s).
            if not isinstance(pump, bool) and not isinstance(turbine, bool):
                self.W_s_turbine, self.W_s_pump = turbine, pump
                self.solve = lambda Q: self._energy_balance(self.W_s_turbine, self.W_s_pump, Q)
                self.unknowns = 1
            elif not isinstance(pump, bool) and isinstance(turbine, bool):
                self.W_s_pump = pump
                self.solve = lambda W_s_turbine, Q: self._energy_balance(W_s_turbine, self.W_s_pump, Q)
                self.unknowns = 2
            elif isinstance(pump, bool) and not isinstance(turbine, bool):
                self.W_s_turbine = turbine
                self.solve = lambda W_s_pump, Q: self._energy_balance(self.W_s_turbine, W_s_pump, Q)
                self.unknowns = 2

    def _energy_balance(self, W_s_turbine, W_s_pump, Q):
        """
        Returns an expression that is equal to zero.

        :param W_s_turbine: The work per unit mass of the turbine(s). Can be a list or float/int.
        :param W_s_pump: The work per unit mass of the pump(s). Can be a list or float/int.
        :param Q: The flow rate. Must be an int/float
        :return: Expression that is equivalent to zero
        """
        self.W_s_turbine, self.W_s_pump, self.Q = W_s_turbine, W_s_pump, Q
        return (self.P_in - self.P_out)/(self.rho*self.g) + (self.V_in**2 - self.V_out**2)/(2*self.g) + self.z_in - \
            self.z_out - np.sum(W_s_turbine)/self.g + np.sum(W_s_pump)/self.g - self.h_f(Q)/self.g

    def _pump_Ws(self, Q):
        """
        Returns the work per unit mass (W_s) for a given flow rate.

        :param Q: The flow rate
        :return: W_s
        """
        return (self.P_out - self.P_in)/self.rho + (self.V_out**2 - self.V_in**2)/2 + self.g*(
                self.z_out - self.z_in) + self.h_f(Q)

    def _turbine_Ws(self, Q):
        """
        Returns the work per unit mass (W_s) for a given flow rate.

        :param Q: The flow rate
        :return: W_s
        """
        return (self.P_in - self.P_out)/self.rho + (self.V_in**2 - self.V_out**2)/2 + self.g*(
                self.z_in - self.z_out) - self.h_f(Q)

    def W_s(self, Q):
        """
        Returns the head of a pump or turbine given the flow rate assuming that the pipe only contains either a pump or
        a turbine. The head is also known as the work per unit mass.

        :param Q: The flow rate
        :return: W_s
        """
        assert self.unknowns == 2, 'You can only acquire the head directly as a function of the flow rate, if there ' \
                                   'is only a pump or only a turbine. If the pipe has both, then the head of either' \
                                   'device must be known (at least one), and solved using the __call__ method and ' \
                                   'fsolve.'

        if self.pump is True:
            return self._pump_Ws(Q)
        else:
            return self._turbine_Ws(Q)

    def delta_P(self, Q):
        """
        Calculates the pressure difference across the pipe at a given flow rate.

        :param Q: The flow rate
        :return: delta_P
        """
        return self.rho*self.W_s(Q)

    def power(self, Q):
        """
        Calculates the power with a given flow rate. If the pipe has a turbine, it calculates the power that gets
        transmitted to the shaft. If the pipe has a pump, then it calculates the power required from the shaft of the
        pump to push the fluid at a given flow rate. This simply means that it will take the efficiency into account.

        :param Q: The flow rate
        :return: power
        """
        power_fluid = self.delta_P(Q)*Q

        if self.turbine:
            return power_fluid*self.eta
        else:
            return power_fluid/self.eta

    def f(self, Q):
        """
        Calculates the friction factor based on the given flow rate. Utilizes the Haaland equations.

        :param Q: The flow rate
        :return: f
        """
        return f(Q, self.D, self.epsilon, self.rho, self.mu)

    def f_T(self):
        """
        Calculates the fully rough friction factor based on the Haaland equation.

        :return: f_T
        """
        return f_T(self.D, self.epsilon)

    def Re(self, Q):
        """
        Calculates Reynold's number. Utilizes the Haaland equation.

        :param Q: The flow rate
        :return: Re
        """
        return Re(Q, self.D, self.rho, self.mu)

    def h_f(self, Q):
        """
        Calculates the losses

        :param Q: The flow rate
        :return: h_f
        """
        return 8*Q**2*(self.C*self.f_T() + self.K + self.L/self.D*self.f(Q))/(np.pi**2*self.D**4)

    def __call__(self, *args):
        """
        This returns an expression that is equal to zero. Used for iterative solutions.

        :param args: Can be 1, 2, or 3 values. If there is only a pump or only a turbine then args is:
                     (W_s_pump, Q) or (W_s_turbine, Q). If there is no pump and no turbine, then args is:
                     (Q, ). If there is a turbine and a pump, then args is:
                     (W_s_turbine, W_s_pump, Q)
        :return: An expression that is equal to zero.
        """
        assert len(args) == self.unknowns, 'Length of args does not match with the number of unknowns.'
        return self.solve(*args)
