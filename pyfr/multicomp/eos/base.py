import numpy as np

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


    @staticmethod
    def poly_reduce(x, y, maxdeg, tol = 5.0):

        ref_poly = np.polynomial.Polynomial.fit(x, y, maxdeg)
        test_x = np.linspace(min(x), max(x), 500)
        ref_y = ref_poly(test_x)
        deg = 0
        while error <= tol or deg <= maxdeg:
            poly = np.polynomial.Polynomial.fit(x, y, deg)
            trial_y = poly(test_x)
            error = np.max(np.abs(trial_y-ref_y)/np.abs(ref_y))

        return poly.convert().coef
