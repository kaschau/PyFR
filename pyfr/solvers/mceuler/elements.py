import numpy as np

from pyfr.solvers.base.elements import ExportableField
from pyfr.solvers.baseadvec import BaseAdvectionElements
from pyfr.multicomp.mcfluid import MCFluidBase, get_mcfluid


class BaseMCFluidElements:
    eos_kernel_module = 'pyfr.solvers.mceuler.kernels.multicomp'

    @classmethod
    def eos_tplargs(cls, ndims, cfg):
        fluid = get_mcfluid(cfg)
        return {
            'ndims': ndims,
            'nvars': len(cls.convars(ndims, cfg)),
            'c': cfg.items_as('constants', float),
            'mcf': fluid,
        }

    @staticmethod
    def privars(ndims, cfg):
        species_names = MCFluidBase.get_species_names(cfg)

        if ndims == 2:
            return ['p', 'u', 'v', 'T'] + species_names[0:-1]
        elif ndims == 3:
            return ['p', 'u', 'v', 'w', 'T'] + species_names[0:-1]

    @staticmethod
    def convars(ndims, cfg):
        species_names = MCFluidBase.get_species_names(cfg)
        if ndims == 2:
            return [f"rho{n}" for n in species_names] + ['rhou', 'rhov', 'E']
        elif ndims == 3:
            return [f"rho{n}" for n in species_names] + ['rhou', 'rhov', 'rhow', 'E']

    dualcoeffs = convars

    @staticmethod
    def visvars(ndims, cfg):
        species_names = MCFluidBase.get_species_names(cfg)
        if ndims == 2:
            varmap = {
                'pressure': ['p'],
                'velocity': ['u', 'v'],
                'temperature': ['T']
            }
        elif ndims == 3:
            varmap = {
                'pressure': ['p'],
                'velocity': ['u', 'v', 'w'],
                'temperature': ['T']
            }
        for sn in species_names[0:-1]:
            varmap[sn] = [sn]

        return varmap

    @staticmethod
    def pri_to_con(pris, cfg):
        fluid = get_mcfluid(cfg)
        return fluid.pri_to_con(pris)

    @staticmethod
    def con_to_pri(cons, cfg):
        fluid = get_mcfluid(cfg)
        return fluid.con_to_pri(cons)

    @staticmethod
    def diff_con_to_pri(cons, diff_cons, cfg):
        fluid = get_mcfluid(cfg)
        return fluid.diff_con_to_pri(cons, diff_cons)

    @staticmethod
    def validate_formulation(ctrl):
        shock_capturing = ctrl.cfg.get('solver', 'shock-capturing', 'none')
        if shock_capturing == 'entropy-filter':
            if ctrl.formulation == 'dual':
                raise ValueError('Entropy filtering not compatible with '
                                 'dual time stepping.')
            elif ctrl.controller_has_variable_dt:
                raise ValueError('Entropy filtering not compatible with '
                                 'adaptive time stepping.')

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Register wavespeed kernel for CFL-based time stepping
        self._be.pointwise.register(
            'pyfr.solvers.mceuler.kernels.wavespeed'
        )

        if self.cfg.getbool('multi-component', 'chemistry', default=False):
            sub_steps = self.cfg.get('multi-component', 'sub-steps', default=0)

            # Sub-stepped chemistry integrates over a compile-time dt and
            # is path dependent; only the instantaneous source is
            # consistent with the implicit residual and its Jacobians
            formulation = self.cfg.get('solver-time-integrator',
                                       'formulation', 'std')
            if formulation == 'implicit' and str(sub_steps) != '0':
                raise ValueError('Implicit time stepping requires '
                                 '[multi-component] sub-steps = 0')

            chem_tplargs = {
                'ndims': self.ndims,
                'nvars': self.nvars,
                'c': self.cfg.items_as('constants', float),
                'mcf': self.mcfluid,
                'dt': self.cfg.getfloat('solver-time-integrator', 'dt'),
            }

            if sub_steps == 'auto':
                max_subs = self.cfg.getfloat('multi-component', 'max-subs', default=10)
                chem_tplargs['max_subs'] = max_subs
                self.add_src_macro('pyfr.solvers.mceuler.kernels.multicomp.chem.finite-rate-auto',
                                   'finite_rate_auto',
                                   chem_tplargs,
                                   False,
                                   True)
            elif not int(sub_steps):
                self.add_src_macro('pyfr.solvers.mceuler.kernels.multicomp.chem.finite-rate',
                                   'finite_rate',
                                   chem_tplargs,
                                   False,
                                   True)
            else:
                chem_tplargs['sub_steps'] = int(sub_steps)
                self.add_src_macro('pyfr.solvers.mceuler.kernels.multicomp.chem.finite-rate-substep',
                                   'finite_rate_substep',
                                   chem_tplargs,
                                   False,
                                   True)

    def init_wavespeed(self):
        self._wspd = self._be.matrix((1, self.neles), tags={'align'})
        self.kernels['wavespeed'] = lambda uin: self._wavespeed_kernel(uin)

        def cfl_getter():
            with np.errstate(divide='ignore'):
                wspd = (2*self.basis.order + 1)*self._wspd.get()[0]
                return np.nan_to_num(1/wspd, posinf=0.0)

        self.export_fields.append(ExportableField(
            name='dt-cfl', shape=(), getter=cfl_getter
        ))
        return self._wspd

    def _wavespeed_kernel(self, uin):
        r, s = self.mesh_regions, self._slice_mat

        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'mcf': self.mcfluid,
            'jac_exprs': self.basis.jac_exprs
        }

        wkerns = []
        for rgn in ('curved', 'linear'):
            if rgn not in r:
                continue

            if rgn == 'curved':
                kw = {'smats': self.curved_smat_at('upts'),
                      'rcpdjac': self.rcpdjac_at('upts', 'curved')}
            else:
                kw = {'verts': self.ploc_at('linspts', 'linear'),
                      'upts': self.upts}

            wkerns.append(self._be.kernel(
                'wavespeed', tplargs=tplargs | {'ktype': rgn},
                dims=[self.nupts, r[rgn]], u=s(self.scal_upts[uin], rgn),
                wspd=self._wspd, **kw
            ))

        if len(wkerns) > 1:
            return self._be.unordered_meta_kernel(wkerns)
        else:
            return wkerns[0]


class MCEulerElements(BaseMCFluidElements, BaseAdvectionElements):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.mcfluid = get_mcfluid(self.cfg)

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.tflux')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'mcf': self.mcfluid,
            'jac_exprs': self.basis.jac_exprs
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
