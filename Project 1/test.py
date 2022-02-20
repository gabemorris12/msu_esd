from msu_esd import Pipe
from scipy.optimize import fsolve
import itertools
import numpy as np
import csv


def network(x, head):
    Q1_, Q2_, Q3_, Q4_, Q5_, Q6_, Q7_, Q8_, Q9_, Q10_, Q11_ = x

    return [
        Q1_ - Q2_ - Q3_,
        Q6_ - Q2_ - Q4_,
        Q3_ - Q4_ - Q5_,
        Q5_ - Q10_ - Q7_,
        Q7_ - Q8_ - Q9_,
        Q11_ - Q6_ - Q7_,
        K_ahu1*Q2_*abs(Q2_) + p2.h(Q2_) + p6.h(Q6_) + p11.h(Q11_) + p1.h(Q1_) - head + 0.1*Q1_*abs(Q1_),
        p4.h(Q4_) + K_ahu3*Q4_*abs(Q4_) - p2.h(Q2_) - K_ahu1*Q2_*abs(Q2_) + p3.h(Q3_),
        p7.h(Q7_) + p8.h(Q8_) + K_ahu4*Q8_*abs(Q8_) - p6.h(Q6_) - p4.h(Q4_) - K_ahu3*Q4_*abs(Q4_) + p5.h(Q5_),
        p9.h(Q9_) + K_ahu4*Q9_*abs(Q9_) - p8.h(Q8_) - K_ahu4*Q8_*abs(Q8_),
        p10.h(Q10_) + K_ahu2*Q10_*abs(Q10_) - p11.h(Q11_) - p9.h(Q9_) - K_ahu4*Q9_*abs(Q9_) - p7.h(Q7_)
    ]


K_ahu1, K_ahu2, K_ahu3, K_ahu4 = 4.5, 4.5, 4.5, 10

rho = 47.3/32.174
mu = 2.72/(3600*32)
epsilon = 0.0005

D = np.array([8.407, 10.482, 12.438, 13.688, 15.67, 17.67, 19.624])/12

p1 = Pipe(D[-1], 2840, epsilon, rho, mu, K=2.58, C=8)
p2 = Pipe(D[1], 2380, epsilon, rho, mu, K=3.93, C=8)
p3 = Pipe(D[2], 1300, epsilon, rho, mu, K=3.18, C=8)
p4 = Pipe(D[1], 1630, epsilon, rho, mu, K=3.33, C=8)
p5 = Pipe(D[2], 3000, epsilon, rho, mu, K=4.68, C=8)
p6 = Pipe(D[2], 5000, epsilon, rho, mu, K=4.68, C=8)
p7 = Pipe(D[2], 1580, epsilon, rho, mu, K=3.38, C=8)
p8 = Pipe(D[0], 1550, epsilon, rho, mu, K=3.33, C=8)
p9 = Pipe(D[0], 1550, epsilon, rho, mu, K=5.28, C=8)
p10 = Pipe(D[1], 5130, epsilon, rho, mu, K=6.78, C=8)
p11 = Pipe(D[2], 1875, epsilon, rho, mu, K=3.33, C=8)

pipes = [eval(f'p{i}') for i in range(1, 12)]

Q_guess = np.array([5, 1, 4, 2, 2, 3, 1, 3, -2, 1, 4])

head_test_file = open('Tests/Overall.csv', 'w', newline='')

for pump_head in range(50, 211, 10):
    print('[Starting]', pump_head)
    products = itertools.product(D, repeat=4)
    tests = []
    count = 0
    for diameters in products:
        count += 1

        d1, d2, d3, d4 = diameters
        p1.D = d1
        p3.D, p5.D, p6.D, p7.D, p11.D = d2, d2, d2, d2, d2
        p2.D, p4.D, p10.D = d3, d3, d3
        p8.D, p9.D = d4, d4

        # noinspection PyTupleAssignmentBalance
        Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11 = fsolve(network, Q_guess, args=(pump_head, ))
        Q = [eval(f'Q{i}') for i in range(1, 12)]

        test = [diameters]
        # The order of the list is:
        # 1) All the velocities are under 9 ft/s
        # 2) 2.67 <= Q2 < 2.77
        # 3) 2.67 <= Q4 < 2.77
        # 4) 1.25 <= Q8 < 1.3
        # 5) 1.25 <= Q9 < 1.3
        # 6) 2.67 <= Q10 < 2.77

        # Test the velocities
        velocity_test = all([pipe.V(Q_value) < 9 for pipe, Q_value in zip(pipes, Q)])
        test.append(velocity_test)

        # Test the flow rates
        Q2_test = 2.67 <= Q2 < 2.77
        Q4_test = 2.67 <= Q4 < 2.77
        Q8_test = 1.25 <= Q8 < 1.3
        Q9_test = 1.25 <= Q9 < 1.3
        Q10_test = 2.67 <= Q10 < 2.77
        for test_ in [Q2_test, Q4_test, Q8_test, Q9_test, Q10_test]:
            test.append(test_)

        tests.append(test)

    with open(f'Tests/Test{pump_head}.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(tests)

    head_tests = [f'{pump_head} ft']
    gt_3, gt_4, gt_5 = [], [], []
    for i, test in enumerate(tests):
        if test.count(True) > 3:
            gt_3.append(i)

        if test.count(True) > 4:
            gt_4.append(i)

        if test.count(True) > 5:
            gt_5.append(i)

    for gt in [gt_3, gt_4, gt_5]:
        head_tests.append(gt)

    if any(gt_3):
        print(f'{len(gt_3)} indices greater than 3 true values.')

    if any(gt_4):
        print(f'{len(gt_4)} indices greater than 4 true values.')

    if any(gt_5):
        print(f'{len(gt_5)} indices greater than 5 true values.')

    head_writer = csv.writer(head_test_file)
    head_writer.writerow(head_tests)

    print('[Ending]', pump_head)

head_test_file.close()
