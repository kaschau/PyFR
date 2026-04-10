from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)
from pyfr.multicomp.mcfluid import get_mcfluid


class TplargsMixin:
    needs_transport = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        if self.cfg.get('solver', 'shock-capturing', 'none') == 'entropy-filter':
            self.d_min = self.cfg.getfloat('solver-entropy-filter', 'd-min',
                                           1e-6)
            self.inte_min = self.cfg.getfloat('solver-entropy-filter', 'inte-min',
                                           1e-6)
        else:
            self.d_min = self.cfg.getfloat('solver-interfaces', 'd-min',
                                           5*self._be.fpdtype_eps)
            self.inte_min = self.cfg.getfloat('solver-interfaces', 'inte-min',
                                           5*self._be.fpdtype_eps)

        self.mcfluid = get_mcfluid(self.cfg, needs_transport=self.needs_transport)

        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             mcf=self.mcfluid,
                             rsolver=rsolver, c=self.c,
                             d_min=self.d_min, inte_min=self.inte_min)

    def validate_species(self):
        sp_names = self.mcfluid.sp_names
        total = sum(float(self.c[n].strip('()')) for n in sp_names[:-1])

        if total > 1.0:
            raise ValueError('BC species mass fractions sum to > 1')
        self.c[sp_names[-1]] = f'({1.0 - total})'


class MCEulerIntInters(TplargsMixin, BaseAdvectionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.intcflux')

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self.scal_lhs, ur=self.scal_rhs, nl=self._pnorm_lhs
        )


class MCEulerMPIInters(TplargsMixin, BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.mpicflux')

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', self._tplargs, dims=[self.ninterfpts],
            ul=self.scal_lhs, ur=self.scal_rhs, nl=self._pnorm_lhs
        )


class MCEulerBaseBCInters(TplargsMixin, BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.bccflux')

        self._tplargs |= dict(bctype=self.type, ninters=self.ninters)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self.scal_lhs, nl=self._pnorm_lhs,
            **self._external_vals
        )

    def comm_entropy_kernel(self, entmin_lhs):
        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.bccent')

        return lambda: self._be.kernel(
            'bccent', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, entmin_lhs=entmin_lhs,
            nl=self._pnorm_lhs, ul=self.scal_lhs, **self._external_vals
        )


class MCEulerSupInflowBCInters(MCEulerBaseBCInters):
    type = 'sup-in-fa'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.mcfluid.sp_names
        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default)
        self.validate_species()


class MCEulerSupOutflowBCInters(MCEulerBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class MCEulerSubOutflowBCInters(MCEulerBaseBCInters):
    type = 'sub-out-fp'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(['p'], lhs)


class MCEulerSlpAdiaWallBCInters(MCEulerBaseBCInters):
    type = 'slp-adia-wall'


class MCEulerConstantMassFlowBCInters(MCEulerBaseBCInters):
    type = 'sub-in-mdot'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.mcfluid.sp_names
        bcvars = ['T', 'mdot-per-area'] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default=default)
        self.validate_species()


class MCEulerCharRiemInvBCInters(MCEulerBaseBCInters):
    type = 'char-riem-inv'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.mcfluid.sp_names
        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default=default)
        self.validate_species()