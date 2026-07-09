from pyfr.fluids import get_fluid
from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.mceuler.elements import BaseMCFluidElements


class MCNavierStokesElements(BaseMCFluidElements,
                             BaseAdvectionDiffusionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        self._add_chem_src()

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        kprefix = 'pyfr.solvers.mcnavstokes.kernels'
        self._be.pointwise.register(f'{kprefix}.tflux')

        fluid = get_fluid(self.cfg, self.ndims)

        if not fluid.provides('mu'):
            raise ValueError('mcnavstokes requires a transport model; set '
                             'transport = constant-props')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'ns': fluid.ns,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'interp_expr': self.basis.interp_expr,
            'fluid': fluid
        }

        # Helpers
        r, s = self.mesh_regions, self._slice_mat

        # Mode-dependent setup
        if self.grad_fusion:
            pts, fused = 'upts', True
            kname = 'tdisf_fused'
        elif 'flux' in self.antialias:
            pts, fused = 'qpts', False
            kname = 'tdisf'
        else:
            pts, fused = 'upts', False
            kname = 'tdisf'

        npts = self.nqpts if pts == 'qpts' else self.nupts
        has_uin = pts == 'upts'

        # Build per-region kernel args
        tdisf = []
        for rgn in ('curved', 'linear'):
            if rgn not in r:
                continue

            ktype = f'{rgn}-fused' if fused else rgn

            # Region-specific geometry kwargs
            kw = {}
            if rgn == 'curved':
                kw['smats'] = self.curved_smat_at(pts)
                if fused:
                    kw['rcpdjac'] = self.rcpdjac_at('upts', 'curved')
            else:
                kw['verts'] = self.ploc_at('linspts', 'linear')

            kw['upts'] = getattr(self, pts)

            if has_uin:
                kw['f'] = s(self._vect_upts, rgn)
                if fused:
                    kw['gradu'] = s(self._grad_upts, rgn)
            else:
                kw['u'] = s(self._scal_qpts, rgn)
                kw['f'] = s(self._vect_qpts, rgn)

            tdisf.append((ktype, r[rgn], rgn, kw))

        if has_uin:
            def tdisf_k(uin):
                return self._make_sliced_kernel(
                    self._be.kernel(
                        'tflux', tplargs=tplargs | {'ktype': kt},
                        dims=[npts, n], u=s(self.scal_upts[uin], rgn), **kw
                    )
                    for kt, n, rgn, kw in tdisf
                )
        else:
            def tdisf_k():
                return self._make_sliced_kernel(
                    self._be.kernel(
                        'tflux', tplargs=tplargs | {'ktype': kt},
                        dims=[npts, n], **kw
                    )
                    for kt, n, rgn, kw in tdisf
                )

        self.kernels[kname] = tdisf_k
