from pyfr.fluids import get_fluid
from pyfr.fluids.constants import RU
from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)


class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        self.fluid = get_fluid(self.cfg, self.ndims)

        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             ns=self.fluid.ns, rsolver=rsolver, c=self.c,
                             fluid=self.fluid)

    def validate_species(self):
        sp_names = self.fluid.sp_names

        try:
            total = sum(float(self.c[n].strip('()')) for n in sp_names[:-1])
        except ValueError:
            exprs = ' - '.join(f'({self.c[n]})' for n in sp_names[:-1])
            self.c[sp_names[-1]] = f'(1.0 - {exprs})'
            return

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
    # Fluids this boundary condition supports as currently implemented
    eos_compat = ('mc-cpg', 'mc-tpg')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.fluid.name not in self.eos_compat:
            raise ValueError(f'Boundary condition {self.type!r} does not '
                             f'support eos {self.fluid.name!r}')

        self._be.pointwise.register('pyfr.solvers.mceuler.kernels.bccflux')

        fluid = self.fluid
        self._tplargs |= dict(
            bctype=self.type, ninters=self.ninters,
            sp_names=fluid.sp_names,
            sp_R=[RU/sp.MW for sp in fluid.species],
            sp_cp0=[getattr(sp, 'cp0', 0.0) for sp in fluid.species]
        )

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self.scal_lhs, nl=self._pnorm_lhs,
            **self._external_vals
        )

    def comm_entropy_kernel(self, entmin_lhs):
        # Physics-specific callback for entropy filtering
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

        sp_names = self.fluid.sp_names
        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default)
        self.validate_species()


class MCEulerSupOutflowBCInters(MCEulerBaseBCInters):
    type = 'sup-out-fn'


class MCEulerSubOutflowBCInters(MCEulerBaseBCInters):
    type = 'sub-out-fp'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(['p'], lhs)


class MCEulerConstantMassFlowBCInters(MCEulerBaseBCInters):
    type = 'sub-in-mdot'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.fluid.sp_names
        bcvars = ['T', 'mdot-per-area'] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default)
        self.validate_species()


class MCEulerCharRiemInvBCInters(MCEulerBaseBCInters):
    type = 'char-riem-inv'

    # Isentropic-relation math assumes a calorically perfect mixture
    eos_compat = ('mc-cpg',)

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.fluid.sp_names
        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default)
        self.validate_species()


class MCEulerSlpAdiaWallBCInters(MCEulerBaseBCInters):
    type = 'slp-adia-wall'
