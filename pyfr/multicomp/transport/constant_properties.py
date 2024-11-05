from pyfr.multicomp.transport.base import BaseTransport


class ConstantProperties(BaseTransport):
    name = 'constant-props'

    def __init__(self, cfg):
        super().__init__(cfg)

        self.input_props = {
            'mu0': None,
            'kappa0': None,
            'Le': None,
        }

        self.consts = {
            'mu0': None,
            'kappa0': None,
            'Le': None,
        }

    @staticmethod
    def compute_consts(props, consts, eos):
        consts['mu0'] = props['mu0']
        consts['kappa0'] = props['kappa0']
        consts['Le'] = props['Le']
