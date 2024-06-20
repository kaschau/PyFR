from pyfr.multicomp.eos.baseEOS import BaseEOS


class cpgEOS(BaseEOS):
    name = 'cpg'
    def __init__(self):
        super().__init__()

        self.input_props = {
            'MW': None,
            'cp0': None,
        }

        self.consts = {
            'MW': None,
            'cp0': None,
        }

    @staticmethod
    def compute_consts(props, consts):
        consts['MW'] = props['MW']
        consts['cp0'] = props['cp0']