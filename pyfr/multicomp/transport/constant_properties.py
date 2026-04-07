from pyfr.multicomp.transport.base import BaseTransport


class ConstantProperties(BaseTransport):
    name = 'constant-props'

    def __init__(self, cfg):
        super().__init__(cfg)

        self.input_props = [
            'mu0',
            'kappa0',
            'Le',
        ]

    def compute_consts(self, props, consts):
        self.consts = consts
        consts['mu0'] = props['mu0']
        consts['kappa0'] = props['kappa0']
        consts['Le'] = props['Le']
