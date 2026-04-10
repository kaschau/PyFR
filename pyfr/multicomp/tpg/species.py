import numpy as np

from pyfr.multicomp import AVOGADRO, EPS0, RU, clean_csigns
from pyfr.multicomp.base import BaseSpecies


class TPGSpecies(BaseSpecies):
    def __init__(self, index, name, spc_sect):
        super().__init__(index, name, spc_sect)

        thermo = spc_sect['thermo']
        self._raw_model = thermo['model'].lower()

        if self._raw_model not in ('nasa7', 'nasa9'):
            raise ValueError(
                f"Species '{name}': unknown thermo model "
                f"'{self._raw_model}'"
            )

        self._raw_ranges = [float(t) for t in thermo['temperature-ranges']]
        self._raw_coeffs = [[float(v) for v in c] for c in thermo['data']]

        # Active thermo data — defaults to raw, EOS class may overwrite
        self.thermo_ranges = list(self._raw_ranges)
        self.thermo_coeffs = [list(c) for c in self._raw_coeffs]

    def bind_strict(self):
        self.eval_mode = 'strict'
        h_idx = 5 if self._raw_model == 'nasa7' else 7
        self.h_ref = self.thermo_coeffs[0][h_idx] * RU / self.MW
        if self._raw_model == 'nasa7':
            self.cp_expr = self._cp_nasa7
            self.h_expr = self._h_nasa7
            self.s_expr = self._s_nasa7
            self.gbs_expr = self._gbs_nasa7
        elif self._raw_model == 'nasa9':
            self.cp_expr = self._cp_nasa9
            self.h_expr = self._h_nasa9
            self.s_expr = self._s_nasa9
            self.gbs_expr = self._gbs_nasa9

    def bind_fast(self):
        self.eval_mode = 'fast'
        self.h_ref = self.thermo_coeffs[0][-2] * RU / self.MW
        self.cp_expr = self._cp_fast
        self.h_expr = self._h_fast
        self.s_expr = self._s_fast
        self.dcp_expr = self._dcp_fast
        self.gbs_expr = self._gbs_fast

    def _multirange(self, var, fn):
        exprs = [fn(var, c) for c in self.thermo_coeffs]
        if len(exprs) == 1:
            return exprs[0]
        result = exprs[-1]
        for i in range(len(exprs) - 2, -1, -1):
            result = f'({var} < {self.thermo_ranges[i+1]}) ? {exprs[i]} : {result}'
        return result

    def _s(self, coeffs, scale=True):
        s = RU / self.MW if scale else 1.0
        return [c * s for c in coeffs]

    def _cp_fast(self, var='T', scale=True):
        c = self.thermo_coeffs[0]
        return self.horner(var, self._s(c[:-2], scale))

    def _h_fast(self, var='T', scale=True):
        c = self.thermo_coeffs[0]
        s = RU / self.MW if scale else 1.0
        return self.horner_integrated(var, self._s(c[:-2], scale),
                                       const=c[-2] * s)

    def _s_fast(self, var='T', scale=True):
        c = self.thermo_coeffs[0]
        sc = self._s(c[:-2], scale)
        s = RU / self.MW if scale else 1.0
        log_term = f'{sc[0]}*log({var})'
        poly = self.horner(var, [sc[i]/i for i in range(1, len(sc))])
        return f'({log_term} + {var}*{poly} + {c[-1] * s})'

    def _dcp_fast(self, var='T', scale=True):
        c = self.thermo_coeffs[0]
        return self.horner_derivative(var, self._s(c[:-2], scale))

    def _cp_nasa7(self, var='T', scale=True):
        return self._multirange(var,
            lambda v, c: self.horner(v, self._s(c[:5], scale)))

    def _h_nasa7(self, var='T', scale=True):
        s = RU / self.MW if scale else 1.0
        return self._multirange(var,
            lambda v, c: self.horner_integrated(v, self._s(c[:5], scale),
                                                const=c[5] * s))

    def _s_nasa7(self, var='T', scale=True):
        def single(v, c):
            sc = self._s(c[:5], scale)
            s = RU / self.MW if scale else 1.0
            log_term = f'{sc[0]}*log({v})'
            poly = self.horner(v, [sc[i]/i for i in range(1, 5)])
            return f'({log_term} + {v}*{poly} + {c[6] * s})'
        return self._multirange(var, single)

    @clean_csigns
    def _gbs_fast(self, var='T', var_logT='logT', var_Tinv='Tinv'):
        c = self.thermo_coeffs[0]
        pc = c[:-2]
        h_const, s_const = c[-2], c[-1]
        terms = [f'{pc[0]} * (1.0 - {var_logT})']
        terms.append(f'{h_const} * {var_Tinv}')
        if len(pc) > 1:
            div_coeffs = [pc[i] / (i * (i + 1)) for i in range(1, len(pc))]
            poly = self.horner(var, [-c for c in div_coeffs])
            terms.append(f'{var}*{poly}')
        terms.append(str(-s_const))
        return '(' + ' + '.join(terms) + ')'

    @clean_csigns
    def _gbs_nasa7(self, var='T', var_logT='logT', var_Tinv='Tinv'):
        def single(v, c):
            pc = c[:5]
            h_const, s_const = c[5], c[6]
            terms = [f'{pc[0]} * (1.0 - {var_logT})']
            terms.append(f'{h_const} * {var_Tinv}')
            div_coeffs = [pc[i] / (i * (i + 1)) for i in range(1, 5)]
            poly = self.horner(v, [-c for c in div_coeffs])
            terms.append(f'{v}*{poly}')
            terms.append(str(-s_const))
            return '(' + ' + '.join(terms) + ')'
        return self._multirange(var, single)

    @clean_csigns
    def _gbs_nasa9(self, var='T', var_logT='logT', var_Tinv='Tinv'):
        def single(v, c):
            terms = [f'-{c[0]}/(2*{v}*{v})']
            terms.append(f'{c[1]}*(1.0 + {var_logT})/{v}')
            terms.append(f'{c[2]} * (1.0 - {var_logT})')
            div_coeffs = [c[i] / ((i - 2) * (i - 1)) for i in range(3, 7)]
            poly = self.horner(v, [-cc for cc in div_coeffs])
            terms.append(f'{v}*{poly}')
            terms.append(f'{c[7]} * {var_Tinv}')
            terms.append(str(-c[8]))
            return '(' + ' + '.join(terms) + ')'
        return self._multirange(var, single)

    @clean_csigns
    def _cp_nasa9(self, var='T', scale=True):
        def single(v, c):
            sc = self._s(c[:7], scale)
            return f'({sc[0]}/({v}*{v}) + {sc[1]}/{v} + {self.horner(v, sc[2:7])})'
        return self._multirange(var, single)

    @clean_csigns
    def _h_nasa9(self, var='T', scale=True):
        s = RU / self.MW if scale else 1.0
        def single(v, c):
            sc = self._s(c[:7], scale)
            neg = f'-{sc[0]}/{v} + {sc[1]}*log({v})'
            pos = self.horner_integrated(v, sc[2:7], const=c[7] * s)
            return f'({neg} + {pos})'
        return self._multirange(var, single)

    @clean_csigns
    def _s_nasa9(self, var='T', scale=True):
        s = RU / self.MW if scale else 1.0
        def single(v, c):
            sc = self._s(c[:7], scale)
            neg = f'-{sc[0]}/(2*{v}*{v}) - {sc[1]}/{v}'
            log_term = f'{sc[2]}*log({v})'
            poly = self.horner(v, [sc[i]/(i-2) for i in range(3, 7)])
            return f'({neg} + {log_term} + {v}*{poly} + {c[8] * s})'
        return self._multirange(var, single)


class TPGTransportSpecies(TPGSpecies):
    def __init__(self, index, name, spc_sect):
        super().__init__(index, name, spc_sect)

        trans = spc_sect['transport']
        self.geometry = trans['geometry']
        self.well_depth = float(trans['well-depth'])
        self.diameter = float(trans['diameter'])
        self.dipole = float(trans.get('dipole', 0.0))
        self.polarizability = float(trans.get('polarizability', 0.0))
        self.rot_relax = float(trans.get('rotational-relaxation', 0.0))
        self.acentric = float(trans.get('acentric-factor', 0.0))

        self.mass = self.MW / AVOGADRO
        self.polar = self.dipole > 0.0

        if self.geometry == 'atom':
            self.rot_dof = 0.0
        elif self.geometry == 'linear':
            self.rot_dof = 1.0
        else:
            self.rot_dof = 1.5

        self.deltastar_self = (
            0.5 * self.dipole**2
            / (4*np.pi * EPS0 * self.well_depth * self.diameter**3)
        )

    def mu_expr(self, var='logT'):
        return self.horner(var, self.muPoly)

    def kappa_expr(self, var='logT'):
        return self.horner(var, self.kappaPoly)

    def dij_expr(self, other_idx, var='logT'):
        return self.horner(var, self.DijPoly[other_idx])
