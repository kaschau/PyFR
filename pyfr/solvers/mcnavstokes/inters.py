import numpy as np

from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)
from pyfr.solvers.mceuler.inters import (MCFluidIntIntersMixin,
                                         MCFluidMPIIntersMixin)
from pyfr.solvers.navstokes.inters import NSCBCMixin
from pyfr.multicomp.mcfluid import MCFluid
from itertools import chain


class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        if shock_capturing == 'entropy-filter':
            self.d_min = self.cfg.getfloat('solver-entropy-filter', 'd-min',
                                           1e-6)
            self.inte_min = self.cfg.getfloat('solver-entropy-filter', 'inte-min',
                                           1e-6)
        else:
            self.d_min = self.cfg.getfloat('solver-interfaces', 'd-min',
                                           5*self._be.fpdtype_eps)
            self.inte_min = self.cfg.getfloat('solver-interfaces', 'inte-min',
                                           5*self._be.fpdtype_eps)
        mcfluid = MCFluid(self.cfg)
        self.c |= mcfluid.consts

        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             rsolver=rsolver,
                             eos = mcfluid.eos, trans = mcfluid.trans,
                             mixing_rule = mcfluid.mixing_rule,
                             shock_capturing=shock_capturing, c=self.c,
                             d_min=self.d_min, inte_min=self.inte_min,
                             visc_sponge=self.visc_sponge)

    def validate_species(self):
        Y = []
        for n in self.c['names']:
            Y.append(float(self.c[n].replace('(','').replace(')','')))

        test = sum(Y)
        if test == 0.0:
            self.c[self.c['names'][-1]] = '(1.)'
        elif test > 1.0:
            raise ValueError('Invalid BC species mass fraction specification.')
        elif test < 1.0:
            self.c[self.c['names'][-1]] = f'({1.0 - test})'

class MCNavierStokesIntInters(TplargsMixin,
                              MCFluidIntIntersMixin,
                              BaseAdvectionDiffusionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.intconu')
        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.intcflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'intconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs,
            ulout=self._comm_lhs, urout=self._comm_rhs
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            nl=self._pnorm_lhs, extrns=self._external_args,
            **self._external_vals
        )


class MCNavierStokesMPIInters(TplargsMixin,
                              MCFluidMPIIntersMixin,
                              BaseAdvectionDiffusionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.mpiconu')
        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.mpicflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'mpiconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs, ulout=self._comm_lhs
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            nl=self._pnorm_lhs, externs=self._external_args,
            **self._external_vals
        )


class MCNavierStokesBaseBCInters(TplargsMixin, BaseAdvectionDiffusionBCInters):
    cflux_state = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Additional BC specific template arguments
        self._tplargs['bctype'] = self.type
        self._tplargs['bccfluxstate'] = self.cflux_state

        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.bcconu')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'bcconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self._scal_lhs,
            ulout=self._comm_lhs, nlin=self._pnorm_lhs,
            **self._external_vals
        )
        if 'nscbc' not in self.type:
            self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.bccflux')
            self.kernels['comm_flux'] = lambda: self._be.kernel(
                'bccflux', tplargs=self._tplargs, dims=[self.ninterfpts],
                extrns=self._external_args, ul=self._scal_lhs,
                gradul=self._vect_lhs, nl=self._pnorm_lhs,
                artviscl=self._artvisc_lhs, **self._external_vals
            )

        if self._ef_enabled:
            self._be.pointwise.register(
                'pyfr.solvers.mcnavstokes.kernels.bccent'
            )

            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'bccent', tplargs=self._tplargs, dims=[self.ninterfpts],
                extrns=self._external_args, entmin_lhs=self._entmin_lhs,
                nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
            )

class MCNavierStokesNoSlpAdiaWallBCInters(MCNavierStokesBaseBCInters):
    type = 'no-slp-adia-wall'
    cflux_state = 'ghost-imperm'


class MCNavierStokesSlpAdiaWallBCInters(MCNavierStokesBaseBCInters):
    type = 'slp-adia-wall'
    cflux_state = None


class MCNavierStokesConstantMassFlowBCInters(MCNavierStokesBaseBCInters):
    type = 'sub-in-mdot'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        bcvars = ['T', 'mdot-per-area']
        bcvars += self.c['names']

        default = {spn: 0 for spn in self.c['names']}

        self.c |= self._exp_opts(bcvars, lhs, default=default)
        self.validate_species()


class MCNavierStokesNoSlpIsotWallBCInters(MCNavierStokesBaseBCInters):
    type = 'no-slp-isot-wall'
    cflux_state = 'ghost-imperm'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(['T'], lhs)
        self.c |= self._exp_opts('uvw'[:self.ndims], lhs,
                                 default={'u': 0, 'v': 0, 'w': 0})


class MCNavierStokesSubOutflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sub-out-fp'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(['p'], lhs)


class MCNavierStokesSupInflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sup-in-fa'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2]
        bcvars += self.c['names']

        default = {spn: 0 for spn in self.c['names']}
        self.c |= self._exp_opts(bcvars, lhs, default=default)
        self.validate_species()


class MCNavierStokesSupOutflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'

class MCNavierStokesCharRiemInvBCInters(MCNavierStokesBaseBCInters):
    type = 'char-riem-inv'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2]
        bcvars += self.c['names']
        default = {spn: 0 for spn in self.c['names']}

        self.c |= self._exp_opts(bcvars, lhs, default=default)
        self.validate_species()

class MCNSCBCSubOutFpBCInters(NSCBCMixin, MCNavierStokesBaseBCInters):

    type = 'sub-out-nscbc-fp'
    decomp_type = 'normal'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        for shape, fidx, lhs_idx in self.ef_pairs:
            # Generate lhs for element-face pair
            lhs_efp = [lhs[i] for i in lhs_idx]
            self.c |= self._exp_opts_ele(['p'], lhs_efp,
                                         self._external_args_efp[shape][fidx],
                                         self._external_vals_efp[shape][fidx])
        self.c['K_p'] = self.cfg.getfloat(cfgsect, 'K_p', default=1.0)

class MCNSCBCSubInFtvyBCInters(NSCBCMixin, MCNavierStokesBaseBCInters):

    type = 'sub-in-nscbc-ftvy'
    decomp_type = 'cartesian'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)
        bcvars = ['T', 'u', 'v', 'w'][:self.ndims + 1]
        for i in chain(bcvars,['Y']):
            self.c[f'K_{i}'] = self.cfg.getfloat(cfgsect, f'K_{i}', default=1.0)

        bcvars += self.c['names']
        default = {spn: 0 for spn in self.c['names']}
        for shape, fidx, lhs_idx in self.ef_pairs:
            # Generate lhs for element-face pair
            lhs_efp = [lhs[i] for i in lhs_idx]

            self.c |= self._exp_opts_ele(bcvars, lhs_efp,
                                         self._external_args_efp[shape][fidx],
                                         self._external_vals_efp[shape][fidx],
                                         default=default)
        self.validate_species()

class MCNSCBCSubInNRIBCInters(NSCBCMixin, MCNavierStokesBaseBCInters):

    type = 'sub-in-nscbc-nri'
    decomp_type = 'normal'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)
        bcvars = ['T', 'un'] + self.c['names']
        force = ['u_a', 'du_a_dt', 'u_v', 'du_v_dt']
        default = {spn: 0 for spn in self.c['names']}
        for shape, fidx, lhs_idx in self.ef_pairs:
            # Generate lhs for element-face pair
            lhs_efp = [lhs[i] for i in lhs_idx]

            self.c |= self._exp_opts_ele(bcvars, lhs_efp,
                                         self._external_args_efp[shape][fidx],
                                         self._external_vals_efp[shape][fidx],
                                         default=default)
            self.c |= self._exp_opts_ele(force, lhs_efp,
                                         self._external_args_efp[shape][fidx],
                                         self._external_vals_efp[shape][fidx],
                                         default={f: 0.0 for f in force})
        for i in ['ac','ut']:
            self.c[f'K_{i}'] = self.cfg.getfloat(cfgsect, f'K_{i}', default=1.0)

        self.validate_species()