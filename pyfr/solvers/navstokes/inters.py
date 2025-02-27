import numpy as np

from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)
from pyfr.solvers.euler.inters import (FluidIntIntersMixin,
                                       FluidMPIIntersMixin)
from pyfr.util import first
from collections import defaultdict


class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        if shock_capturing == 'entropy-filter':
            self.p_min = self.cfg.getfloat('solver-entropy-filter', 'p-min',
                                           1e-6)
        else:
            self.p_min = self.cfg.getfloat('solver-interfaces', 'p-min',
                                           5*self._be.fpdtype_eps)

        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             rsolver=rsolver, visc_corr=visc_corr,
                             shock_capturing=shock_capturing, c=self.c,
                             p_min=self.p_min)


class NavierStokesIntInters(TplargsMixin,
                            FluidIntIntersMixin,
                            BaseAdvectionDiffusionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.intconu')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.intcflux')

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


class NavierStokesMPIInters(TplargsMixin,
                            FluidMPIIntersMixin,
                            BaseAdvectionDiffusionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.mpiconu')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.mpicflux')

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


class NavierStokesBaseBCInters(TplargsMixin, BaseAdvectionDiffusionBCInters):
    cflux_state = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Additional BC specific template arguments
        self._tplargs['bctype'] = self.type
        self._tplargs['bccfluxstate'] = self.cflux_state

        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.bcconu')
        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.bccflux')

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
                'pyfr.solvers.navstokes.kernels.bccent'
            )

            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'bccent', tplargs=self._tplargs, dims=[self.ninterfpts],
                extrns=self._external_args, entmin_lhs=self._entmin_lhs,
                nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
            )


class NavierStokesNoSlpIsotWallBCInters(NavierStokesBaseBCInters):
    type = 'no-slp-isot-wall'
    cflux_state = 'ghost-imperm'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c['cpTw'], = self._eval_opts(['cpTw'])
        self.c |= self._exp_opts('uvw'[:self.ndims], lhs,
                                 default={'u': 0, 'v': 0, 'w': 0})


class NavierStokesNoSlpAdiaWallBCInters(NavierStokesBaseBCInters):
    type = 'no-slp-adia-wall'
    cflux_state = 'ghost-imperm'


class NavierStokesSlpAdiaWallBCInters(NavierStokesBaseBCInters):
    type = 'slp-adia-wall'
    cflux_state = None


class NavierStokesCharRiemInvBCInters(NavierStokesBaseBCInters):
    type = 'char-riem-inv'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )


class NavierStokesSupInflowBCInters(NavierStokesBaseBCInters):
    type = 'sup-in-fa'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )


class NavierStokesSupOutflowBCInters(NavierStokesBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class NavierStokesSubInflowFrvBCInters(NavierStokesBaseBCInters):
    type = 'sub-in-frv'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(
            ['rho', 'u', 'v', 'w'][:self.ndims + 1], lhs,
            default={'u': 0, 'v': 0, 'w': 0}
        )


class NavierStokesSubInflowFtpttangBCInters(NavierStokesBaseBCInters):
    type = 'sub-in-ftpttang'
    cflux_state = 'ghost'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        gamma = self.cfg.getfloat('constants', 'gamma')

        # Pass boundary constants to the backend
        self.c['cpTt'], = self._eval_opts(['cpTt'])
        self.c['pt'], = self._eval_opts(['pt'])
        self.c['Rdcp'] = (gamma - 1.0)/gamma

        # Calculate u, v velocity components from the inflow angle
        theta = self._eval_opts(['theta'])[0]*np.pi/180.0
        velcomps = np.array([np.cos(theta), np.sin(theta), 1.0])

        # Adjust u, v and calculate w velocity components for 3-D
        if self.ndims == 3:
            phi = self._eval_opts(['phi'])[0]*np.pi/180.0
            velcomps[:2] *= np.sin(phi)
            velcomps[2] *= np.cos(phi)

        self.c['vc'] = velcomps[:self.ndims]


class NavierStokesSubOutflowBCInters(NavierStokesBaseBCInters):
    type = 'sub-out-fp'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['p'], lhs)


class NavierStokesCharacteristicBoundaryCondition(NavierStokesBaseBCInters):

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.bccflux_nscbc')

        self._tplargs_efp = defaultdict(dict)

        self._scal_upts = defaultdict(dict)
        self._scal_fpts = defaultdict(dict)
        self._vect_upts = defaultdict(dict)
        self._vect_fpts = defaultdict(dict)
        self._pnorm_facefpts = defaultdict(dict)
        self._smats_upts = defaultdict(dict)
        self._jacs_facefpts = defaultdict(dict)

        # lhs length
        self._dim_lhs = defaultdict(dict)

        self.kernels['comm_flux'] = lambda: self.gen_nscbc_kerns()

        # Create required element-face pairs
        ef_pairs = []
        for shape in set(t[0] for t in lhs):
            for fidx in range(len(self.elemap[shape].basis.faces)):
                lhs_efp = [t for t in lhs if t[0] == shape and t[2] == fidx]
                if lhs_efp:
                    ef_pairs.append((shape, fidx, lhs_efp))


        for shape, fidx, lhs_efp in ef_pairs:
            basis = self.elemap[shape].basis
            nupts = basis.nupts
            nfpts = basis.nfpts
            nfacefpts = basis.nfacefpts[fidx]
            ndims = self.ndims
            nvars = self.nvars

            # Generate lhs for element-face pair
            self._dim_lhs[shape][fidx] = len(lhs_efp)

            # Store tplargs for this element-face pair
            self._tplargs_efp[shape][fidx] = self._tplargs.copy()
            tplargs_efp = self._tplargs_efp[shape][fidx]

            facefpts = basis.facefpts[fidx]
            tplargs_efp['m0'] = basis.m0[facefpts]
            tplargs_efp['m1'] = basis.m1.reshape(nupts,ndims,nupts)
            tplargs_efp['m2'] = basis.m2.reshape(nfpts,ndims,nupts)[facefpts]
            tplargs_efp['m11'] = basis.m11[facefpts, facefpts]
            tplargs_efp['m12'] = basis.m12[facefpts]

            tplargs_efp['nupts'] = nupts
            tplargs_efp['nfpts'] = nfpts
            tplargs_efp['nfacefpts'] = nfacefpts
            tplargs_efp['facefpts'] = basis.facefpts[fidx]
            tplargs_efp['bnorm_facefpts'] = basis.norm_fpts[basis.facefpts[fidx]]

            scal_upts = self._scal_upts_view(lhs_efp, '_get_scal_upts_for_inter_ele')
            self._scal_upts[shape][fidx] = scal_upts

            scal_fpts = self._scal_fpts_view(lhs_efp, '_get_scal_fpts_for_inter_ele')
            self._scal_fpts[shape][fidx] = scal_fpts

            vect_upts = self._vect_upts_view(lhs_efp, '_get_vect_upts_for_inter_ele')
            self._vect_upts[shape][fidx] = vect_upts

            vect_fpts = self._vect_fpts_view(lhs_efp, '_get_vect_fpts_for_inter_ele')
            self._vect_fpts[shape][fidx] = vect_fpts

            pnorm_facefpts = self._fwise_const_mat(lhs_efp, '_get_pnorms_facefpts')
            self._pnorm_facefpts[shape][fidx] = pnorm_facefpts

            smats_upts = self._ewise_const_mat(lhs_efp, '_get_smats_upts')
            self._smats_upts[shape][fidx] = smats_upts

            jacs_facefpts = self._fwise_const_mat(lhs_efp, '_get_jacs_facefpts')
            self._jacs_facefpts[shape][fidx] = jacs_facefpts

    def gen_nscbc_kerns(self):
        kerns = []
        for shape in self._tplargs_efp.keys():
            for fidx in self._tplargs_efp[shape].keys():

                tplargs_efp = self._tplargs_efp[shape][fidx]

                kerns.append(self._be.kernel(
                    'bccflux_nscbc', tplargs=tplargs_efp,
                    dims=[self._dim_lhs[shape][fidx]],
                    extrns=self._external_args,
                    u_upts=self._scal_upts[shape][fidx],
                    u_fpts=self._scal_fpts[shape][fidx],
                    gradu_upts=self._vect_upts[shape][fidx],
                    gradu_fpts=self._vect_fpts[shape][fidx],
                    nl_ffpt=self._pnorm_facefpts[shape][fidx],
                    smats_upts=self._smats_upts[shape][fidx],
                    jacs_ffpt=self._jacs_facefpts[shape][fidx],
                    **self._external_vals))

        return self._be.unordered_meta_kernel(kerns)


class NSCBCSubOutFPInters(NavierStokesCharacteristicBoundaryCondition):

    type = 'sub-out-nscbc-fp'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['p'], lhs)
        self.c['sigma'] = self.cfg.getfloat(cfgsect, 'sigma', default=0.25)
