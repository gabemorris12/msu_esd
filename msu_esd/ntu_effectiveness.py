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


def parallel_single_pass(x, C, find='NTU'):
    """
    For a parallel flow single pass

    :param x: NTU or effectiveness. If find="NTU", then x is expected to be the effectiveness. x can be an array.
    :param C: C_min/C_max
    :param find: Specify whether to find the NTU or effectiveness
    :return: NTU or effectiveness
    """
    eff_func = lambda ntu: (1 - np.exp(-ntu*(C + 1)))/(C + 1)
    ntu_func = lambda e: np.log(-1/(e*(C + 1) - 1))/(C + 1)
    return _finder([ntu_func, eff_func], x, find=find)


def counter_single_pass(x, C, find='NTU'):
    """
    For a counter flow single pass

    :param x: NTU or effectiveness. If find="NTU", then x is expected to be the effectiveness. x can be an array.
    :param C: C_min/C_max
    :param find: Specify whether to find the NTU or effectiveness
    :return: NTU or effectiveness
    """
    if C == 1 and find in ['NTU', 'ntu']:
        return x/(1 - x)
    elif C == 1 and find in ['e', 'epsilon', 'effectiveness']:
        return x/(1 + x)
    elif C == 1:
        raise Exception('Valid arguments for "find" are "NTU" or "effectiveness", "e", and "epsilon."')
    eff_func = lambda ntu: (1 - np.exp(-ntu*(1 - C)))/(-C*np.exp(-ntu*(1 - C)) + 1)
    ntu_func = lambda e: -np.log((C*e - 1)/(e - 1))/(C - 1)
    return _finder([ntu_func, eff_func], x, find=find)


def shell_and_tube_one_shell_pass(x, C, find='NTU'):
    """
    For a shell and tube with one shell pass and some multiple of two tube passes

    :param x: NTU or effectiveness. If find="NTU", then x is expected to be the effectiveness. x can be an array.
    :param C: C_min/C_max
    :param find: Specify whether to find the NTU or effectiveness
    :return: NTU or effectiveness
    """
    eff_func = lambda ntu: 2/(
            C + 1 + (1 + np.exp(-ntu*(C**2 + 1)**0.5))*(C**2 + 1)**0.5/(1 - np.exp(-ntu*(C**2 + 1)**0.5)))
    ntu_func = lambda e: np.log((C*e - e*np.sqrt(C**2 + 1.0) + e - 2.0)/(C*e + e*(C**2 + 1.0)**0.5 + e - 2.0))/np.sqrt(
        C**2 + 1.0)
    return _finder([ntu_func, eff_func], x, find=find)


def shell_and_tube_n_shell_passes():
    """
    For a shell and tube with n number of shell passes and some multiple 2n two tube passes
    """
    raise NotImplemented('Function not yet implemented')


def cross_flow_unmixed(x, C, find='NTU'):
    """
    For a cross flow with both streams unmixed

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


def cross_flow_both_mixed():
    raise NotImplemented('Function not yet implemented')


def cross_flow_C_min_unmixed():
    raise NotImplemented('Function not yet implemented')


def cross_flow_C_max_unmixed():
    raise NotImplemented('Function not yet implemented')
