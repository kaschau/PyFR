import numpy as np

class BaseEOS:
    name = None

    def __init__(self):
        pass

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
