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

        assert not all([pump, turbine]), 'As of right now, pipe objects cannot have a pump and turbine.'
        assert any([pump, turbine]), 'As of right now, pipe objects must have a pump or turbine.'

        if pump:
            self.Ws = self._pump_Ws
        else:
            self.Ws = self._turbine_Ws

    def _pump_Ws(self, Q):
        """
        Returns the work per unit mass (Ws) for a given flow rate.

        :param Q: The flow rate
        :return: Ws
        """
        return (self.P_out - self.P_in)/self.rho + (self.V_out**2 - self.V_in**2)/2 + self.g*(
                self.z_out - self.z_in) + self.h_f(Q)

    def _turbine_Ws(self, Q):
        """
        Returns the work per unit mass (Ws) for a given flow rate.

        :param Q: The flow rate
        :return: Ws
        """
        return (self.P_in - self.P_out)/self.rho + (self.V_in**2 - self.V_out**2)/2 + self.g*(
                self.z_in - self.z_out) - self.h_f(Q)

    def delta_P(self, Q):
        """
        Calculates the pressure difference across the pipe at a given flow rate.

        :param Q: The flow rate
        :return: delta_P
        """
        return self.rho*self.Ws(Q)

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

    def find_flow_or_work(self, Ws, Q):
        """
        Returns an expression that is equal to zero. This is used for solving for the flow rate or the work per unit
        mass of the fluid.

        :param Ws: Work per unit mass of the fluid
        :param Q: The flow rate
        :return: Ws - Ws(Q), an expression that is equal to zero.
        """
        return Ws - self.Ws(Q)

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

    def __call__(self, Ws, Q):
        return self.find_flow_or_work(Ws, Q)
