from lzma import MF_BT2
import numpy as np

from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)
from pyfr.solvers.mceuler.inters import (MCFluidIntIntersMixin,
                                         MCFluidMPIIntersMixin)
from pyfr.multicomp.mcfluid import MCFluid


class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        if shock_capturing == 'entropy-filter':
            self.Y_tol = self.cfg.getfloat('solver-entropy-filter', 'Y-tol',
                                           1e-12)
            self.p_min = self.cfg.getfloat('solver-entropy-filter', 'p-min',
                                           1e-6)
        else:
            self.Y_tol = self.cfg.getfloat('solver-interfaces', 'Y-tol',
                                           5*self._be.fpdtype_eps)
            self.p_min = self.cfg.getfloat('solver-interfaces', 'p-min',
                                           5*self._be.fpdtype_eps)
        mcfluid = MCFluid(self.cfg)
        self.c |= mcfluid.consts

        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             rsolver=rsolver, visc_corr=visc_corr,
                             eos = mcfluid.eos, trans = mcfluid.trans,
                             shock_capturing=shock_capturing, c=self.c,
                             Y_tol=self.Y_tol, p_min=self.p_min)


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
            nl=self._pnorm_lhs
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
            nl=self._pnorm_lhs
        )


class MCNavierStokesBaseBCInters(TplargsMixin, BaseAdvectionDiffusionBCInters):
    cflux_state = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Additional BC specific template arguments
        self._tplargs['bctype'] = self.type
        self._tplargs['bccfluxstate'] = self.cflux_state

        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.bcconu')
        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.bccflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'bcconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self._scal_lhs,
            ulout=self._comm_lhs, nlin=self._pnorm_lhs,
            **self._external_vals
        )
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
            self._tplargs['e_func'] = self.cfg.get('solver-entropy-filter',
                                                   'e-func', 'numerical')

            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'bccent', tplargs=self._tplargs, dims=[self.ninterfpts],
                extrns=self._external_args, entmin_lhs=self._entmin_lhs,
                nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
            )

class MCNavierStokesNoSlpAdiaWallBCInters(MCNavierStokesBaseBCInters):
    type = 'no-slp-adia-wall'
    cflux_state = 'ghost'


class MCNavierStokesSlpAdiaWallBCInters(MCNavierStokesBaseBCInters):
    type = 'slp-adia-wall'
    cflux_state = None


class MCNavierStokesConstantMassFlowBCInters(MCNavierStokesBaseBCInters):
    type = 'sub-in-mdot'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        bcvars = ['T', 'mdot-per-area']
        bcvars += self.c['names'][:-1]

        default = {spn: 0 for spn in self.c['names'][:-1]}

        self.c |= self._exp_opts(bcvars, lhs, default=default)


class MCNavierStokesNoSlpIsotWallBCInters(MCNavierStokesBaseBCInters):
    type = 'no-slp-isot-wall'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['T'], lhs)
        self.c |= self._exp_opts('uvw'[:self.ndims], lhs,
                                 default={'u': 0, 'v': 0, 'w': 0})


class MCNavierStokesSubOutflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sub-out-fp'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['p'], lhs)


class MCNavierStokesSupInflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sup-in-fa'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2]
        bcvars += self.c['names'][:-1]

        default = {spn: 0 for spn in self.c['names'][:-1]}
        self.c |= self._exp_opts(bcvars, lhs, default=default)


class MCNavierStokesSupOutflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'
