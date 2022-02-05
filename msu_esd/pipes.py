
class Pipe:
    def __init__(self, D, L, epsilon, rho, mu, K=0, C=0, P_in=0, P_out=0, V_in=0, V_out=0, z_in=0, z_out=0):
        self.D, self.L, self.epsilon, self.rho, self.mu = D, L, epsilon, rho, mu
        self.K, self.C = K, C
        self.P_in, self.P_out = P_in, P_out
        self.V_in, self.V_out = V_in, V_out
        self.z_in, self.z_out = z_in, z_out
