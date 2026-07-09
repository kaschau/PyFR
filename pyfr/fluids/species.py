import functools
import re

from pyfr.fluids.constants import RU


# Tidy up +- sign sequences in generated C expressions
def clean_csigns(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        s = fn(*args, **kwargs)
        s = re.sub(r'\+\s*\-', '- ', s)
        s = re.sub(r'\-\s*\-', '+ ', s)
        s = re.sub(r'\-\s*\+', '- ', s)
        s = re.sub(r'\+\s*\+', '+ ', s)
        s = re.sub(r'=\s*\+\s*', '= ', s)
        return s
    return wrapper


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

        return ('(' + str(coeffs[0]) + ' + ' + var + '*'
                + BaseSpecies.horner(var, coeffs[1:]) + ')')

    @staticmethod
    def horner_integrated(var, coeffs, const=0):
        int_coeffs = [c / (i + 1) for i, c in enumerate(coeffs)]
        expr = f'{var}*{BaseSpecies.horner(var, int_coeffs)}'
        if const:
            expr = f'({expr} + {const})'

        return expr


class CPGSpecies(BaseSpecies):
    def __init__(self, index, name, spc_sect):
        super().__init__(index, name, spc_sect)
        self.cp0 = float(spc_sect['cp0'])


class TPGSpecies(BaseSpecies):
    def __init__(self, index, name, spc_sect):
        super().__init__(index, name, spc_sect)

        thermo = spc_sect['thermo']
        if thermo['model'].lower() != 'nasa7':
            raise ValueError(f'Species {name!r}: only NASA7 thermo is '
                             'currently supported')

        self.thermo_ranges = [float(t) for t in thermo['temperature-ranges']]
        self.thermo_coeffs = [[float(v) for v in c] for c in thermo['data']]

    def _scaled(self, coeffs):
        return [c*RU/self.MW for c in coeffs]

    def _multirange(self, var, fn):
        exprs = [fn(var, c) for c in self.thermo_coeffs]
        if len(exprs) == 1:
            return exprs[0]

        result = exprs[-1]
        for i in range(len(exprs) - 2, -1, -1):
            result = (f'({var} < {self.thermo_ranges[i + 1]})'
                      f' ? {exprs[i]} : {result}')

        return result

    # Per-mass specific heat and enthalpy as C expressions in var
    def cp_expr(self, var='T'):
        return self._multirange(
            var, lambda v, c: self.horner(v, self._scaled(c[:5]))
        )

    def h_expr(self, var='T'):
        return self._multirange(
            var, lambda v, c: self.horner_integrated(
                v, self._scaled(c[:5]), const=c[5]*RU/self.MW
            )
        )

    def s_expr(self, var='T'):
        def single(v, c):
            sc = self._scaled(c)

            return (f'({sc[0]}*log({v})'
                    f' + {v}*{self.horner(v, [sc[i]/i for i in range(1, 5)])}'
                    f' + {sc[6]})')

        return self._multirange(var, single)

    # Non-dimensional Gibbs free energy G/(R T) for equilibrium constants
    @clean_csigns
    def gbs_expr(self, var='T', var_logT='logT', var_Tinv='Tinv'):
        def single(v, c):
            pc = c[:5]
            h_const, s_const = c[5], c[6]
            terms = [f'{pc[0]} * (1.0 - {var_logT})']
            terms.append(f'{h_const} * {var_Tinv}')
            div_coeffs = [pc[i] / (i * (i + 1)) for i in range(1, 5)]
            poly = self.horner(v, [-cc for cc in div_coeffs])
            terms.append(f'{v}*{poly}')
            terms.append(str(-s_const))

            return '(' + ' + '.join(terms) + ')'

        return self._multirange(var, single)
