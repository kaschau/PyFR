import numpy as np

class BaseEOS:
    name = None

    def __init__(self):
        pass

    def validate_Y_ics(self, Yk):
        for Y in Yk:
            if np.max(Y) > 1.0:
                raise ValueError('Species mass > 1.0 detected in ICs')
            elif np.min(Y) < 0.0:
                raise ValueError('Species mass < 0.0 detected in ICs')
