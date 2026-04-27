import numpy as np

from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)
from pyfr.solvers.mceuler.inters import TplargsMixin as _MCEulerTplargsMixin
from pyfr.solvers.navstokes.inters import NSCBCMixin


class MCNSCBCMixin(NSCBCMixin):
    _nscbc_kern = 'pyfr.solvers.mcnavstokes.kernels.bccflux_nscbc'


class TplargsMixin(_MCEulerTplargsMixin):
    needs_transport = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._tplargs['shock_capturing'] = self.cfg.get(
            'solver', 'shock-capturing', 'none'
        )

class MCNavierStokesIntInters(TplargsMixin, BaseAdvectionDiffusionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.intconu')
        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.intcflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'intconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self.scal_lhs, urin=self.scal_rhs,
            ulout=self._comm_lhs, urout=self._comm_rhs
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self.scal_lhs, ur=self.scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artvisc=self.artvisc, nl=self._pnorm_lhs
        )


class MCNavierStokesMPIInters(TplargsMixin, BaseAdvectionDiffusionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.mpiconu')
        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.mpicflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'mpiconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self.scal_lhs, urin=self.scal_rhs, ulout=self._comm_lhs
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self.scal_lhs, ur=self.scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artvisc=self.artvisc, nl=self._pnorm_lhs
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
            extrns=self._external_args, ulin=self.scal_lhs,
            ulout=self._comm_lhs, nlin=self._pnorm_lhs,
            **self._external_vals
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self.scal_lhs,
            gradul=self._vect_lhs, nl=self._pnorm_lhs,
            artvisc=self.artvisc, **self._external_vals
        )

    def comm_entropy_kernel(self, entmin_lhs):
        self._be.pointwise.register('pyfr.solvers.mcnavstokes.kernels.bccent')

        return lambda: self._be.kernel(
            'bccent', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, entmin_lhs=entmin_lhs,
            nl=self._pnorm_lhs, ul=self.scal_lhs, **self._external_vals
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

        sp_names = self.mcfluid.sp_names
        bcvars = ['T', 'mdot-per-area'] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

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

        sp_names = self.mcfluid.sp_names
        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

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

        sp_names = self.mcfluid.sp_names
        bcvars = ['T', 'p', 'u', 'v', 'w'][:self.ndims + 2] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        self.c |= self._exp_opts(bcvars, lhs, default=default)
        self.validate_species()


class MCNSCBCSubOutFpBCInters(MCNSCBCMixin, MCNavierStokesBaseBCInters):

    type = 'sub-out-nscbc-fp'
    decomp_type = 'normal'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        for etype, fidx in self.ef_pairs:
            lhs_efp = self._lhs_efp[etype][fidx]
            self.c |= self._exp_opts_ele(
                ['p'], lhs_efp,
                self._external_args_efp[etype][fidx],
                self._external_vals_efp[etype][fidx],
            )
        self.c['K_p'] = self.cfg.getfloat(cfgsect, 'K_p', default=0.25)


class MCNSCBCSubInFtvyBCInters(MCNSCBCMixin, MCNavierStokesBaseBCInters):

    type = 'sub-in-nscbc-ftvy'
    decomp_type = 'cartesian'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.mcfluid.sp_names
        bcvars = ['T', 'u', 'v', 'w'][:self.ndims + 1] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        for etype, fidx in self.ef_pairs:
            lhs_efp = self._lhs_efp[etype][fidx]
            self.c |= self._exp_opts_ele(
                bcvars, lhs_efp,
                self._external_args_efp[etype][fidx],
                self._external_vals_efp[etype][fidx],
                default=default,
            )

        for i in ['T', 'u', 'v', 'w'][:self.ndims + 1]:
            self.c[f'K_{i}'] = self.cfg.getfloat(cfgsect, f'K_{i}', default=0.25)
        self.c['K_Y'] = self.cfg.getfloat(cfgsect, 'K_Y', default=0.25)
        self.validate_species()


class MCNSCBCSubInNRIBCInters(MCNSCBCMixin, MCNavierStokesBaseBCInters):

    type = 'sub-in-nscbc-nri'
    decomp_type = 'normal'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        sp_names = self.mcfluid.sp_names
        force = ['u_a', 'du_a_dt', 'u_v', 'du_v_dt']

        bcvars = ['T', 'un'] + list(sp_names)
        default = {spn: 0 for spn in sp_names}

        for etype, fidx in self.ef_pairs:
            lhs_efp = self._lhs_efp[etype][fidx]
            self.c |= self._exp_opts_ele(
                bcvars, lhs_efp,
                self._external_args_efp[etype][fidx],
                self._external_vals_efp[etype][fidx],
                default=default,
            )
            self.c |= self._exp_opts_ele(
                force, lhs_efp,
                self._external_args_efp[etype][fidx],
                self._external_vals_efp[etype][fidx],
                default={f: 0.0 for f in force},
            )

        for i in ['ac', 'ut']:
            self.c[f'K_{i}'] = self.cfg.getfloat(cfgsect, f'K_{i}', default=0.25)
        self.validate_species()