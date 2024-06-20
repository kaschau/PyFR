from pyfr.multicomp.eos.baseEOS import BaseEOS


class tpgEOS(BaseEOS):
    name = 'tpg'
    def __init__(self):
        super().__init__()

        self.input_props = {
            'MW': None,
            'NASA7': None,
        }

        self.consts = {
            'MW': None,
            'NASA7': None,
        }

    @staticmethod
    def compute_consts(props, consts):
        consts['MW'] = props['MW']
        consts['NASA7'] = props['NASA7']