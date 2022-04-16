from msu_esd import Pipe
import numpy as np
import csv
from scipy.optimize import fsolve

# Define known constants
K1, K2, K3, K4 = 4.5, 4.5, 4.5, 10  # Loss coefficients

rho = 48.05/32.174  # In slugs per cubic feet
mu = 3.715/(3600*32.174)  # In slugs per (ft s) or lbf*s per ft squared
epsilon = 0.0005  # In ft

D20, D12, D10, D8 = np.array([18.812, 11.938, 10.02, 7.981])/12  # Diameters in ft

# Possible main line diameters in ft
D_main = [1.567667, 1.885333, 2.552000, 2.718667, 2.875000, 3.375000]
# Possible supply line diameters in ft
D_supply = [0.665083, 0.835000, 0.994833, 1.093667, 1.250000, 1.406333, 1.567667]
# Possible ahu line diameters in ft for line 2, 8, 9, 10
D_air_handling = [0.420583, 0.505417, 0.665083, 0.835000, 0.994833, 1.093667, 1.250000, 1.406333]
# Possible ahu line diameters in ft for line 4
D_line_4 = [0.665083, 0.835000, 0.994833, 1.093667, 1.250000, 1.406333]

p1 = Pipe(D20, 2840, epsilon, rho, mu, C=8)
p2 = Pipe(D10, 2380, epsilon, rho, mu, C=8)
p3 = Pipe(D12, 1300, epsilon, rho, mu, C=8)
p4 = Pipe(D10, 1630, epsilon, rho, mu, C=8)
p5 = Pipe(D12, 3000, epsilon, rho, mu, C=8)
p6 = Pipe(D12, 5000, epsilon, rho, mu, C=8)
p7 = Pipe(D12, 1580, epsilon, rho, mu, C=8)
p8 = Pipe(D8, 1550, epsilon, rho, mu, C=8)
p9 = Pipe(D8, 1550, epsilon, rho, mu, C=8)
p10 = Pipe(D10, 5130, epsilon, rho, mu, C=8)
p11 = Pipe(D12, 1875, epsilon, rho, mu, C=8)

Q_guess = np.array([5, 1, 4, 2, 2, 3, 1, 3, -2, 1, 4])
pipes = [eval(f'p{i}') for i in range(1, 12)]


def unbalanced(x, Ws):
    Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11 = x
    # all expressions need to be set to zero
    return [
        Q1 - Q2 - Q3,
        Q6 - Q2 - Q4,
        Q3 - Q4 - Q5,
        Q5 - Q10 - Q7,
        Q7 - Q8 - Q9,
        Q11 - Q6 - Q7,
        K1*Q2*abs(Q2) + p2.h(Q2) + p6.h(Q6) + p11.h(Q11) + p1.h(Q1) - Ws + 0.1*Q1*abs(Q1),
        p4.h(Q4) + K3*Q4*abs(Q4) - p2.h(Q2) - K1*Q2*abs(Q2) + p3.h(Q3),
        p7.h(Q7) + p8.h(Q8) + K4*Q8*abs(Q8) - p6.h(Q6) - p4.h(Q4) - K3*Q4*abs(Q4) + p5.h(Q5),
        p9.h(Q9) + K4*Q9*abs(Q9) - p8.h(Q8) - K4*Q8*abs(Q8),
        p10.h(Q10) + K2*Q10*abs(Q10) - p11.h(Q11) - p9.h(Q9) - K4*Q9*abs(Q9) - p7.h(Q7)
    ]


def test(Q, diameters):
    Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11 = Q
    velocity_test = all([p.V(Q_) < 9 for p, Q_ in zip(pipes, Q)])
    tests = [diameters, velocity_test]
    Q2_test = 2.713 <= Q2 < 2.808
    Q4_test = 2.713 <= Q4 < 2.808
    Q8_test = 1.272 <= Q8 < 1.316
    Q9_test = 1.272 <= Q9 < 1.316
    Q10_test = 2.713 <= Q10 < 2.808
    [tests.append(t) for t in [Q2_test, Q4_test, Q8_test, Q9_test, Q10_test]]
    return tests


def iteration(Ws):
    """
    Results will display:
    1) Diameters in the following order [Line 1; Line 3, 5, 6, 7, 11; Line 2, 10; Line 8, 9; Line 4]
    2) True for all velocities being under 9
    3) 2.713 <= Q2 < 2.808
    4) 2.713 <= Q4 < 2.808
    5) 1.272 <= Q8 < 1.316
    6) 1.272 <= Q9 < 1.316
    7) 2.713 <= Q10 < 2.808
    """
    file = open(f'Tests/Test{Ws}.csv', 'w', newline='')
    writer = csv.writer(file)
    print(f'Starting {Ws}')
    for d_main in D_main:
        p1.D = d_main
        for d_supply in D_supply:
            p3.D, p5.D, p6.D, p7.D, p11.D = d_supply, d_supply, d_supply, d_supply, d_supply
            for d_ahu in D_air_handling:
                p2.D, p10.D = d_ahu, d_ahu
                for d_line_8_9 in D_air_handling:
                    p8.D, p9.D = d_line_8_9, d_line_8_9
                    for d_line_4 in D_line_4:
                        p4.D = d_line_4
                        sol = fsolve(unbalanced, Q_guess, args=(Ws, ))
                        result = test(sol, [d_main, d_supply, d_ahu, d_line_8_9, d_line_4])
                        writer.writerow(result)
                        if result[-5:].count(True) >= 4:
                            print(Ws, result)
    print(f'Ending {Ws}')

    file.close()


if __name__ == '__main__':
    for i in range(180, 231):
        iteration(i)
