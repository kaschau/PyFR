import time

from pyfr.cache import memoize
from pyfr.integrators.implicit.precond import Preconditioner


class MCChemPreconditioner(Preconditioner):
    """Point-implicit chemistry preconditioner.

    Approximates (I - gdt*df/du)^-1 with the block-diagonal-by-point
    inverse of I - gdt*d(src)/du built from the analytic chemistry
    Jacobian.  Blocks are nvars x nvars per solution point (memory and
    build cost are nupts^2 times smaller than element block-Jacobi) and
    are rebuilt from the current Newton iterate at every construct call,
    which is one pointwise kernel launch.  Captures the chemical
    stiffness only; convective/viscous coupling is left to the Krylov
    solver.
    """

    name = 'point-chem'

    def __init__(self, backend, system, ext, *, fd_eps=0, pcdtype=None):
        super().__init__(backend, system, ext, fd_eps=fd_eps,
                         pcdtype=pcdtype)
        self.computed = False

        self._eles = eles = list(system.ele_map.values())

        cfg = eles[0].cfg
        if not cfg.getbool('multi-component', 'chemistry', False):
            raise ValueError('point-chem preconditioner requires '
                             '[multi-component] chemistry = true')

        backend.pointwise.register(
            'pyfr.solvers.mceuler.kernels.ptchemprecond'
        )
        backend.pointwise.register(
            'pyfr.solvers.mceuler.kernels.ptchemapply'
        )

        # Per-point inverse blocks for each element type
        self._minvs = [
            backend.matrix((e.nupts, e.nvars*e.nvars, e.neles),
                           tags={'align'})
            for e in eles
        ]

    def _tplargs(self, e):
        return {'nupts': e.nupts, 'nvars': e.nvars, 'ndims': e.ndims,
                'mcf': e.mcfluid, 'c': e.cfg.items_as('constants', float)}

    @memoize
    def _get_build_kerns(self, u_reg):
        kerns = []
        for etidx, e in enumerate(self._eles):
            kerns.append(self.backend.kernel(
                'ptchemprecond', tplargs=self._tplargs(e), dims=[e.neles],
                u=self.system.ele_banks[etidx][u_reg],
                minv=self._minvs[etidx]
            ))
        return kerns

    def construct(self, t, u_reg, gamma_dt, rhs_fn, f0_reg, up_reg, add_fn,
                  ptmp, eps_scales=()):
        # The build is a single cheap kernel launch, so always rebuild
        # from the current Newton iterate
        t0 = time.perf_counter()

        kerns = self._get_build_kerns(u_reg)
        for k in kerns:
            k.bind(gdt=gamma_dt)
        self.backend.run_kernels(kerns)

        self.computed = True
        self.gdt_built = gamma_dt
        self.build_t = t
        self.build_wtime = time.perf_counter() - t0
        self.build_wtime_total += self.build_wtime
        self.nbuilds += 1

    def invalidate(self):
        self.computed = False

    def apply_kernel(self, emats, etidx, in_reg, out_reg, in_scale=(),
                     out_scale=()):
        e = self._eles[etidx]
        tplargs = self._tplargs(e) | {'in_scale': in_scale,
                                      'out_scale': out_scale}

        return self.backend.kernel(
            'ptchemapply', tplargs=tplargs, dims=[e.neles],
            x=emats[in_reg], minv=self._minvs[etidx], y=emats[out_reg]
        )
