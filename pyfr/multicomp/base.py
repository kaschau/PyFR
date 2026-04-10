class BaseEos:
    name = None


class BaseSpecies:
    def __init__(self, index, name, spc_sect):
        self.index = int(index)
        self.name = name
        self.composition = dict(spc_sect['composition'])
        self.MW = float(spc_sect['MW'])

    @staticmethod
    def horner(var, coeffs):
        if not coeffs:
            return '0'
        if len(coeffs) == 1:
            return str(coeffs[0])
        return '(' + str(coeffs[0]) + ' + ' + var + '*' + BaseSpecies.horner(var, coeffs[1:]) + ')'

    @staticmethod
    def horner_integrated(var, coeffs, const=0):
        # integral of c0 + c1*x + c2*x^2 + ... = c0*x + c1*x^2/2 + ...
        # = x*(c0 + x*(c1/2 + x*(c2/3 + ...))) + const
        int_coeffs = [c/(i+1) for i, c in enumerate(coeffs)]
        inner = BaseSpecies.horner(var, int_coeffs)
        expr = f'{var}*{inner}'
        if const:
            expr = f'({expr} + {const})'
        return expr

    @staticmethod
    def horner_derivative(var, coeffs):
        dcoeffs = [(i+1)*c for i, c in enumerate(coeffs[1:])]
        return BaseSpecies.horner(var, dcoeffs) if dcoeffs else '0'
