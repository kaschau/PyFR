from functools import cache

import numpy as np

from pyfr.fluids import get_fluid
from pyfr.fluids.constants import RU
from pyfr.fluids.readers.cantera import read_cantera_yaml
from pyfr.solvers.baseadvec import BaseAdvectionElements


@cache
def _species_count(path):
    species_sects, _ = read_cantera_yaml(path)

    return len(species_sects)


class BaseMCFluidElements:
    @classmethod
    def _fluid(cls, cfg, nvals):
        # Both privar and convar counts equal ns + ndims + 1
        path = cfg.getpath('multi-component', 'species')
        ndims = nvals - _species_count(path) - 1

        return get_fluid(cfg, ndims)

    @staticmethod
    def privars(ndims, cfg):
        return get_fluid(cfg, ndims).privars

    @staticmethod
    def convars(ndims, cfg):
        return get_fluid(cfg, ndims).convars

    dualcoeffs = convars

    @staticmethod
    def visvars(ndims, cfg):
        return get_fluid(cfg, ndims).visvars

    @classmethod
    def pri_to_con(cls, pris, cfg):
        return cls._fluid(cfg, len(pris)).pri_to_con(pris)

    @classmethod
    def con_to_pri(cls, cons, cfg):
        return cls._fluid(cfg, len(cons)).con_to_pri(cons)

    def _add_chem_src(self):
        cfg = self.cfg
        if not cfg.getbool('multi-component', 'chemistry', False):
            return

        fluid = get_fluid(cfg, self.ndims)
        prec = cfg.get('backend', 'precision', 'double')
        fpd = np.float32 if prec == 'single' else np.float64

        chem_tplargs = {
            'ns': fluid.ns, 'fluid': fluid, 'RU': RU,
            'dt': cfg.getfloat('solver-time-integrator', 'dt'),
            'participates': fluid.species_participates,
            'fpdtype_min': float(np.finfo(fpd).tiny),
            'fpdtype_eps': float(np.finfo(fpd).eps)
        }

        kpre = 'pyfr.solvers.mceuler.kernels.chem'
        sub_steps = cfg.get('multi-component', 'sub-steps', '0')
        if sub_steps == 'auto':
            chem_tplargs['max_subs'] = cfg.getint('multi-component',
                                                  'max-subs', 10)
            self.add_src_macro(f'{kpre}.finite-rate-auto',
                               'finite_rate_auto', chem_tplargs,
                               False, True)
        elif not int(sub_steps):
            self.add_src_macro(f'{kpre}.finite-rate', 'finite_rate',
                               chem_tplargs, False, True)
        else:
            chem_tplargs['sub_steps'] = int(sub_steps)
            self.add_src_macro(f'{kpre}.finite-rate-substep',
                               'finite_rate_substep', chem_tplargs,
                               False, True)


class MCEulerElements(BaseMCFluidElements, BaseAdvectionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        self._add_chem_src()

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.tflux')

        fluid = get_fluid(self.cfg, self.ndims)

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'ns': fluid.ns,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'fluid': fluid
        }

        # Helpers
        tdisf = []
        c, l = 'curved', 'linear'
        r, s = self.mesh_regions, self._slice_mat
        slicedk = self._make_sliced_kernel

        if c in r and 'flux' not in self.antialias:
            tdisf.append(lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'curved'},
                dims=[self.nupts, r[c]], u=s(self.scal_upts[uin], c),
                f=s(self._vect_upts, c), smats=self.curved_smat_at('upts')
            ))
        elif c in r:
            tdisf.append(lambda: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'curved'},
                dims=[self.nqpts, r[c]], u=s(self._scal_qpts, c),
                f=s(self._vect_qpts, c), smats=self.curved_smat_at('qpts')
            ))

        if l in r and 'flux' not in self.antialias:
            tdisf.append(lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'linear'},
                dims=[self.nupts, r[l]], u=s(self.scal_upts[uin], l),
                f=s(self._vect_upts, l), verts=self.ploc_at('linspts', l),
                upts=self.upts
            ))
        elif l in r:
            tdisf.append(lambda: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'linear'},
                dims=[self.nqpts, r[l]], u=s(self._scal_qpts, l),
                f=s(self._vect_qpts, l), verts=self.ploc_at('linspts', l),
                upts=self.qpts
            ))

        if 'flux' not in self.antialias:
            self.kernels['tdisf'] = lambda uin: slicedk(k(uin) for k in tdisf)
        else:
            self.kernels['tdisf'] = lambda: slicedk(k() for k in tdisf)
