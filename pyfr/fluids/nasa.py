import numpy as np


def eval_nasa7_cp(T, a):
    return a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*a[4])))


def eval_nasa7_h(T, a):
    return T*(a[0] + T*(a[1]/2 + T*(a[2]/3 + T*(a[3]/4 + T*a[4]/5)))) + a[5]


def eval_nasa7_s(T, a):
    return (a[0]*np.log(T)
            + T*(a[1] + T*(a[2]/2 + T*(a[3]/3 + T*a[4]/4))) + a[6])


def eval_over_ranges(T, ranges, coeffs, eval_fn):
    T = np.atleast_1d(np.asarray(T, dtype=float))
    result = np.empty_like(T)

    for j in range(len(coeffs)):
        if j == 0:
            mask = T <= ranges[j + 1]
        elif j == len(coeffs) - 1:
            mask = T > ranges[j]
        else:
            mask = (T > ranges[j]) & (T <= ranges[j + 1])

        if np.any(mask):
            result[mask] = eval_fn(T[mask], coeffs[j])

    return result
