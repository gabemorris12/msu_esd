import numpy as np
from scipy.optimize import fsolve


def _finder(funcs, *args, find='NTU'):
    ntu_func, eff_func = funcs

    if find == 'NTU' or find == 'ntu':
        return ntu_func(*args)
    elif find == 'effectiveness' or find == 'e' or find == 'epsilon':
        return eff_func(*args)
    else:
        raise Exception('Valid arguments for "find" are "NTU" or "effectiveness", "e", and "epsilon."')


def cross_flow_unmixed(x, C, find='NTU'):
    """
    Calculate either NTU or effectiveness

    :param x: NTU or effectiveness. If find="NTU", then x is expected to be the effectiveness. x can be an array.
    :param C: C_min/C_max
    :param find: Specify whether to find the NTU or effectiveness
    :return: NTU or effectiveness
    """
    eff_func = lambda ntu: 1 - np.exp(ntu**0.22*(-1 + np.exp(-C*ntu**0.78))/C)

    if find == 'NTU' or find == 'ntu':
        try:
            iter(x)
        except TypeError:
            ntu_func = lambda ntu: x - eff_func(ntu)
            return fsolve(ntu_func, np.array([1, ]))[0]
        else:
            def ntu_func(e):
                ntu_func__ = lambda ntu: e - eff_func(ntu)
                return fsolve(ntu_func__, np.array([0, ]))[0]
            return np.array(list(map(ntu_func, x)))
    elif find == 'effectiveness' or find == 'e' or find == 'epsilon':
        return eff_func(x)
    else:
        raise Exception('Valid arguments for "find" are "NTU" or "effectiveness", "e", and "epsilon."')
