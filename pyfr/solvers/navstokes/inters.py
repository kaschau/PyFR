import numpy as np

from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)
from pyfr.solvers.euler.inters import (FluidIntIntersMixin,
                                       FluidMPIIntersMixin,
                                       MassFlowBCMixin)
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
        self.kernels['con_u'] = lambda: self._be.kernel(
            'bcconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self._scal_lhs,
            ulout=self._comm_lhs, nlin=self._pnorm_lhs,
            **self._external_vals
        )
        if 'nscbc' not in self.type:
            self._be.pointwise.register('pyfr.solvers.navstokes.kernels.bccflux')
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

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

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

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )


class NavierStokesSupInflowBCInters(NavierStokesBaseBCInters):
    type = 'sup-in-fa'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )


class NavierStokesSupOutflowBCInters(NavierStokesBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class NavierStokesSubInflowFrvBCInters(NavierStokesBaseBCInters):
    type = 'sub-in-frv'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

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

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(['p'], lhs)


class NavierStokesCharRiemInvMassFlowBCInters(MassFlowBCMixin,
                                              NavierStokesBaseBCInters):
    type = 'char-riem-inv-mass-flow'
    cflux_state = 'ghost'


class NavierStokesCharacteristicBoundaryCondition(NavierStokesBaseBCInters):

    @staticmethod
    def _newCS_2d(n):
        return np.array([-n[1],n[0]])

    @staticmethod
    def _newCS_3d(n):
        nz = n[2]
        if abs(nz) < 0.9:
            temp = np.array([0, 0, 1])
        else:
            temp = np.array([1, 0, 0])

        t1 = np.cross(temp,n)
        t1 /= np.linalg.norm(t1)

        t2 = np.cross(n, t1)

        return t1, t2

    def newCS(self, n):
        if len(n) == 2:
            return self._newCS_2d(n)
        else:
            return self._newCS_3d(n)

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        kname = 'pyfr.solvers.navstokes.kernels.bccflux_nscbc'
        self._be.pointwise.register(kname)

        self._tplargs_efp = defaultdict(dict)
        self._external_args_efp = defaultdict(dict)
        self._external_vals_efp = defaultdict(dict)

        self._scal_upts = defaultdict(dict)
        self._scal_fpts = defaultdict(dict)
        self._grad_upts = defaultdict(dict)
        self._vect_fpts = defaultdict(dict)
        self._normnl_facefpts = defaultdict(dict)
        self._smats_upts = defaultdict(dict)
        self._jacs_facefpts = defaultdict(dict)

        # lhs length
        self._dim_lhs = defaultdict(dict)

        # Register NSCBC kernel under a different name so it can be scheduled after other BCs
        self.kernels['nscbc_flux'] = lambda: self.gen_nscbc_kerns()

        # Create required element-face pairs
        self.ef_pairs = []
        for shape in set(t[0] for t in lhs):
            basis = self.elemap[shape].basis
            for fidx in range(len(basis.faces)):
                lhs_idx = [i for i,t in enumerate(lhs)
                           if t[0] == shape and t[2] == fidx]
                if lhs_idx:
                    self.ef_pairs.append((shape, fidx, lhs_idx))


        for shape, fidx, lhs_idx in self.ef_pairs:
            self._dim_lhs[shape][fidx] = len(lhs_idx)

            # Generate lhs for element-face pair
            lhs_efp = [lhs[i] for i in lhs_idx]

            # Store tplargs for this element-face pair
            self._tplargs_efp[shape][fidx] = self._tplargs.copy()
            self._external_args_efp[shape][fidx] = self._external_args.copy()
            self._external_vals_efp[shape][fidx] = self._external_vals.copy()

            # Basis in regular orientation
            ele = self.elemap[shape]
            basis = ele.basis
            nupts = basis.nupts
            nfpts = basis.nfpts
            nfacefpts = basis.nfacefpts[fidx]
            ndims = self.ndims
            facefpts = basis.facefpts[fidx]
            nintfpts = nfpts - nfacefpts
            intfpts = [i for i in range(nfpts) if i not in facefpts]

            tplargs_efp = self._tplargs_efp[shape][fidx]

            tplargs_efp['nupts'] = nupts
            tplargs_efp['nfpts'] = nfpts
            tplargs_efp['nfacefpts'] = nfacefpts
            tplargs_efp['nintfacefpts'] = nintfpts
            tplargs_efp['facefpts'] = facefpts
            tplargs_efp['intfpts'] = intfpts
            norms = basis.norm_fpts[facefpts]
            tplargs_efp['magnl'] = np.linalg.norm(norms, axis=-1)
            norms = norms/tplargs_efp['magnl'][:, None]
            if self.normal == 'inward':
                norms *= -1.0
            tplargs_efp['bnorms'] = norms

            tplargs_efp['t1s'] = np.empty(norms.shape)
            if ndims == 3:
                tplargs_efp['t2s'] = np.empty(norms.shape)
            for i, norm in enumerate(norms):
                if ndims == 2:
                    tplargs_efp['t1s'][i] = self.newCS(norm)
                else:
                    tplargs_efp['t1s'][i], tplargs_efp['t2s'][i] = self.newCS(norm)

            tplargs_efp['m2'] = basis.m2.reshape(nfpts,ndims,nupts)
            tplargs_efp['m12'] = basis.m12[facefpts]

            # Compute inverse of correction function matrix G for face flux points
            # GB_ij represents correction function j evaluated at flux point i
            # m11[i,j] is the correction function from flux point j at flux point i
            GB = np.zeros((nfacefpts, nfacefpts))
            for i, fpt_i in enumerate(facefpts):
                for j, fpt_j in enumerate(facefpts):
                    GB[i, j] = basis.m11[fpt_i, fpt_j]
            GI = np.zeros((nfacefpts, nintfpts))
            for i, fpt_i in enumerate(facefpts):
                for j, fpt_j in enumerate(intfpts):
                    GI[i, j] = basis.m11[fpt_i, fpt_j]

            # Compute inverse of G matrix
            GB_inv = np.linalg.inv(GB)
            tplargs_efp['GB_inv'] = GB_inv
            tplargs_efp['GB_inv_GI'] = GB_inv @ GI

            method = '_get_scal_upts_for_inter_ele'
            scal_upts = self._scal_upts_view(lhs_efp, method)
            self._scal_upts[shape][fidx] = scal_upts

            method = '_get_scal_fpts_for_inter_ele'
            scal_fpts = self._scal_fpts_view(lhs_efp, method)
            self._scal_fpts[shape][fidx] = scal_fpts

            method = '_get_grad_upts_for_inter_ele'
            grad_upts = self._grad_upts_view(lhs_efp, method)
            self._grad_upts[shape][fidx] = grad_upts

            method = '_get_vect_fpts_for_inter_ele'
            vect_fpts = self._vect_fpts_view(lhs_efp, method)
            self._vect_fpts[shape][fidx] = vect_fpts

            ## Get normals at all flux points (not just boundary face)
            ## We need a custom method to get normals for all flux points in element
            if self.normal == 'inward':
                method = '_get_inward_normnls_facefpts'
            else:
                method = '_get_normnls_facefpts'
            normnl_facefpts = self._ewise_const_mat(lhs_efp, method)
            self._normnl_facefpts[shape][fidx] = normnl_facefpts

            method = '_get_smats_upts'
            smats_upts = self._ewise_const_mat(lhs_efp, method)
            self._smats_upts[shape][fidx] = smats_upts

            method = '_get_jacs_facefpts'
            jacs_facefpts = self._ewise_const_mat(lhs_efp, method)
            self._jacs_facefpts[shape][fidx] = jacs_facefpts

    def gen_nscbc_kerns(self):
        kerns = []
        for shape in self._tplargs_efp.keys():
            for fidx in self._tplargs_efp[shape].keys():

                kerns.append(self._be.kernel(
                    'bccflux_nscbc',
                    tplargs=self._tplargs_efp[shape][fidx],
                    dims=[self._dim_lhs[shape][fidx]],
                    extrns=self._external_args_efp[shape][fidx],
                    u_upts=self._scal_upts[shape][fidx],
                    u_fpts=self._scal_fpts[shape][fidx],
                    gradu_upts=self._grad_upts[shape][fidx],
                    gradu_fpts=self._vect_fpts[shape][fidx],
                    normnl_ffpts=self._normnl_facefpts[shape][fidx],
                    smats_upts=self._smats_upts[shape][fidx],
                    jacs_ffpts=self._jacs_facefpts[shape][fidx],
                    **self._external_vals_efp[shape][fidx]))

        return self._be.unordered_meta_kernel(kerns)


class NSCBCSubOutFpBCInters(NavierStokesCharacteristicBoundaryCondition):

    type = 'sub-out-nscbc-fp'
    normal = 'outward'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        for shape, fidx, lhs_idx in self.ef_pairs:
            # Generate lhs for element-face pair
            lhs_efp = [lhs[i] for i in lhs_idx]
            self.c |= self._exp_opts_ele(['p'], lhs_efp,
                                         self._external_args_efp[shape][fidx],
                                         self._external_vals_efp[shape][fidx])
        self.c['K_p'] = self.cfg.getfloat(cfgsect, 'K_p', default=1.0)

class NSCBCSubInFrvBCInters(NavierStokesCharacteristicBoundaryCondition):

    type = 'sub-in-nscbc-frv'
    normal = 'inward'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)
        for shape, fidx, lhs_idx in self.ef_pairs:
            # Generate lhs for element-face pair
            lhs_efp = [lhs[i] for i in lhs_idx]

            self.c |= self._exp_opts_ele(
                ['rho', 'u', 'v', 'w'][:self.ndims + 1], lhs_efp,
                self._external_args_efp[shape][fidx],
                self._external_vals_efp[shape][fidx],
            )
        for i in ['rho', 'u', 'v', 'w'][:self.ndims + 1]:
            self.c[f'K_{i}'] = self.cfg.getfloat(cfgsect, f'K_{i}', default=1.0)