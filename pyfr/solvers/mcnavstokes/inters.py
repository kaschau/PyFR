from pyfr.fluids.constants import RU
from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)
from pyfr.solvers.mceuler.inters import TplargsMixin as MCEulerTplargsMixin


class TplargsMixin(MCEulerTplargsMixin):
    pass


class MCNavierStokesIntInters(TplargsMixin, BaseAdvectionDiffusionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        kprefix = 'pyfr.solvers.mcnavstokes.kernels'
        self._be.pointwise.register(f'{kprefix}.intconu')
        self._be.pointwise.register(f'{kprefix}.intcflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'intconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self.scal_lhs, urin=self.scal_rhs,
            ulout=self._comm_lhs, urout=self._comm_rhs
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self.scal_lhs, ur=self.scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs, nl=self._pnorm_lhs
        )


class MCNavierStokesMPIInters(TplargsMixin, BaseAdvectionDiffusionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        kprefix = 'pyfr.solvers.mcnavstokes.kernels'
        self._be.pointwise.register(f'{kprefix}.mpiconu')
        self._be.pointwise.register(f'{kprefix}.mpicflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'mpiconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self.scal_lhs, urin=self.scal_rhs, ulout=self._comm_lhs
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self.scal_lhs, ur=self.scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs, nl=self._pnorm_lhs
        )


class MCNavierStokesBaseBCInters(TplargsMixin, BaseAdvectionDiffusionBCInters):
    cflux_state = None

    # Fluids this boundary condition supports as currently implemented
    eos_compat = ('mc-cpg', 'mc-tpg')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.fluid.name not in self.eos_compat:
            raise ValueError(f'Boundary condition {self.type!r} does not '
                             f'support eos {self.fluid.name!r}')

        kprefix = 'pyfr.solvers.mcnavstokes.kernels'
        self._be.pointwise.register(f'{kprefix}.bcconu')
        self._be.pointwise.register(f'{kprefix}.bccflux')

        fluid = self.fluid
        self._tplargs |= dict(
            bctype=self.type, bccfluxstate=self.cflux_state,
            ninters=self.ninters, sp_names=fluid.sp_names,
            sp_R=[RU/sp.MW for sp in fluid.species],
            sp_cp0=[getattr(sp, 'cp0', 0.0) for sp in fluid.species]
        )

        self.kernels['con_u'] = lambda: self._be.kernel(
            'bcconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self.scal_lhs,
            ulout=self._comm_lhs, nlin=self._pnorm_lhs,
            **self._external_vals
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self.scal_lhs,
            gradul=self._vect_lhs, nl=self._pnorm_lhs,
            **self._external_vals
        )

    def comm_entropy_kernel(self, entmin_lhs):
        # Physics-specific callback for entropy filtering
        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.bccent')

        return lambda: self._be.kernel(
            'bccent', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, entmin_lhs=entmin_lhs,
            nl=self._pnorm_lhs, ul=self.scal_lhs, **self._external_vals
        )


class MCNavierStokesNoSlpAdiaWallBCInters(MCNavierStokesBaseBCInters):
    type = 'no-slp-adia-wall'
    cflux_state = 'ghost-imperm'


class MCNavierStokesNoSlpIsotWallBCInters(MCNavierStokesBaseBCInters):
    type = 'no-slp-isot-wall'
    cflux_state = 'ghost-imperm'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(['T'], lhs)
        self.c |= self._exp_opts('uvw'[:self.ndims], lhs,
                                 default={'u': 0, 'v': 0, 'w': 0})


class MCNavierStokesSlpAdiaWallBCInters(MCNavierStokesBaseBCInters):
    type = 'slp-adia-wall'
    cflux_state = None


class MCNavierStokesSupInflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sup-in-fa'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.fluid.sp_names
        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default)
        self.validate_species()


class MCNavierStokesSupOutflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class MCNavierStokesSubOutflowBCInters(MCNavierStokesBaseBCInters):
    type = 'sub-out-fp'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(['p'], lhs)


class MCNavierStokesConstantMassFlowBCInters(MCNavierStokesBaseBCInters):
    type = 'sub-in-mdot'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.fluid.sp_names
        bcvars = ['T', 'mdot-per-area'] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default)
        self.validate_species()


class MCNavierStokesCharRiemInvBCInters(MCNavierStokesBaseBCInters):
    type = 'char-riem-inv'
    cflux_state = 'ghost'

    # Isentropic-relation math assumes a calorically perfect mixture
    eos_compat = ('mc-cpg',)

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.fluid.sp_names
        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default)
        self.validate_species()
