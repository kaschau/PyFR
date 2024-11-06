import numpy as np

def poly_reduce(x, y, maxdeg, tol = 0.01, w = None):

    ref_poly = np.polynomial.Polynomial.fit(x, y, maxdeg, w=w)
    test_x = np.linspace(min(x), max(x), 500)
    ref_y = ref_poly(test_x)
    deg = 0
    error = np.inf
    while deg <= maxdeg:
        poly = np.polynomial.Polynomial.fit(x, y, deg, w=w)
        trial_y = poly(test_x)
        error = np.max(np.abs(trial_y-ref_y)/np.abs(ref_y))
        deg += 1
        if error <= tol:
            break

    return poly

class BaseEOS:
    name = None

    def __init__(self, cfg):
        self.cfg = cfg

    def validate_Y_ics(self, Yk):
        Ysum = 0.0
        for Y in Yk:
            if np.max(Y) > 1.0:
                raise ValueError('Species mass fraction > 1.0 detected in ICs')
            elif np.min(Y) < 0.0:
                raise ValueError('Species mass fraction < 0.0 detected in ICs')
            Ysum += Y
        if np.max(Ysum) > 1.0:
            raise ValueError('Species mass fraction sum > 1.0 detected in ICs')