from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)
from pyfr.multicomp.mcfluid import MCFluid


class MCFluidIntIntersMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.cfg.get('solver', 'shock-capturing') == 'entropy-filter':
            self._be.pointwise.register('pyfr.solvers.mceuler.kernels.intcent')

            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'intcent', tplargs={}, dims=[self.ninters],
                entmin_lhs=self._entmin_lhs, entmin_rhs=self._entmin_rhs
            )


class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        if self.cfg.get('solver', 'shock-capturing') == 'entropy-filter':
            self.p_min = self.cfg.getfloat('solver-entropy-filter', 'p-min',
                                           1e-6)
        else:
            self.p_min = self.cfg.getfloat('solver-interfaces', 'p-min',
                                           5*self._be.fpdtype_eps)
        self.mcfluid = MCFluid(self.cfg)
        self.c |= self.mcfluid.consts

        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             eos=self.mcfluid.eos,
                             rsolver=rsolver, c=self.c, p_min=self.p_min)


class MCFluidMPIIntersMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.cfg.get('solver', 'shock-capturing') == 'entropy-filter':
            self._be.pointwise.register('pyfr.solvers.mceuler.kernels.mpicent')

            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'mpicent', tplargs={}, dims=[self.ninters],
                entmin_lhs=self._entmin_lhs, entmin_rhs=self._entmin_rhs
            )


class MCEulerIntInters(TplargsMixin, MCFluidIntIntersMixin,
                     BaseAdvectionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.intcflux')

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )


class MCEulerMPIInters(TplargsMixin, MCFluidMPIIntersMixin,
                     BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.mpicflux')

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )


class MCEulerBaseBCInters(TplargsMixin, BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.bccflux')

        self._tplargs |= dict(bctype=self.type, ninters=self.ninters)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs, nl=self._pnorm_lhs,
            **self._external_vals
        )

        if self.cfg.get('solver', 'shock-capturing') == 'entropy-filter':
            self._be.pointwise.register('pyfr.solvers.mceuler.kernels.bccent')

            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'bccent', tplargs=self._tplargs, dims=[self.ninterfpts],
                extrns=self._external_args, entmin_lhs=self._entmin_lhs,
                nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
            )


class MCEulerSupInflowBCInters(MCEulerBaseBCInters):
    type = 'sup-in-fa'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2]
        bcvars += self.c['names'][:-1]

        default = {spn: 0 for spn in self.c['names'][:-1]}

        self.c |= self._exp_opts(bcvars, lhs, default)


class MCEulerSupOutflowBCInters(MCEulerBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class MCEulerSubOutflowBCInters(MCEulerBaseBCInters):
    type = 'sub-out-fp'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['p'], lhs)


class MCEulerSlpAdiaWallBCInters(MCEulerBaseBCInters):
    type = 'slp-adia-wall'


class MCEulerConstantMassFlowBCInters(MCEulerBaseBCInters):
    type = 'sub-in-mdot'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        bcvars = ['T', 'mdot-per-area']
        bcvars += self.c['names'][:-1]

        default = {spn: 0 for spn in self.c['names'][:-1]}

        self.c |= self._exp_opts(bcvars, lhs, default=default)