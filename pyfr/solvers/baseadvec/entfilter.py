import numpy as np

from pyfr.quadrules import get_quadrule
from pyfr.solvers.base.elements import ExportableField


# Reference-space parametric measure of each face kind.  Used to build the
# divergence-theorem integration weight for shock-normal detectors.
_FACE_REF_MEASURE = {
    'line': 2.0,    # s in [-1, 1]
    'tri':  2.0,    # PyFR standard tri area
    'quad': 4.0,    # s, t in [-1, 1]
}

# Reference-space volume of each element type.
_CELL_REF_VOLUME = {
    'quad': 4.0,        # [-1, 1]^2
    'hex':  8.0,        # [-1, 1]^3
    'tri':  2.0,
    'tet':  4.0/3.0,
    'pri':  4.0,        # tri area * line length = 2 * 2
    'pyr':  8.0/3.0,
}


class EntropyFilter:
    def __init__(self, backend, cfg, system, int_inters, mpi_inters,
                 bc_inters):
        self._be = backend

        # Register pointwise kernel templates
        kprefix = f'pyfr.solvers.{system.ef_solver}.kernels'
        backend.pointwise.register(f'{kprefix}.entropylocal')
        backend.pointwise.register(f'{kprefix}.entropyfilter')

        # Per-element-type setup
        self._entmin = {}
        for etype, eles in system.ele_map.items():
            if eles.basis.order > 0:
                self._setup_etype(etype, eles, cfg, system.nonce)

        backend.commit()

        # Create interface views and register comm_entropy kernels
        self._setup_interfaces(system, int_inters, mpi_inters, bc_inters)

    def _setup_etype(self, etype, eles, cfg, nonce):
        be = self._be
        nfaces = len(eles.nfacefpts)
        neles = eles.neles

        # Allocate one minimum entropy value per face
        entmin_int = np.full((nfaces, neles), -be.fpdtype_max,
                             dtype=be.fpdtype)
        entmin = be.matrix((nfaces, neles), tags={'align'},
                            extent=nonce + 'entmin_int',
                            initval=entmin_int)
        self._entmin[etype] = entmin

        # Allocate space for filter strength (1 = no filter, 0 = max)
        ef_filter = be.matrix((1, neles),
                               extent=nonce + 'ef_filter',
                               tags={'align'})

        # Register exportable field for filter strength
        eles.export_fields.append(ExportableField(
            name='ef-filter', shape=(),
            getter=lambda: ef_filter.get()[0]
        ))

        # Per-cell shock-normal (reference-space) and its magnitude.  Written
        # inside the kernel when a cell hits admissibility failure; zero in
        # smooth cells.  Exported for visualisation/diagnostics.
        shock_normal = be.matrix((eles.ndims, neles),
                                  extent=nonce + 'shock_normal',
                                  tags={'align'})
        shock_normal_mag = be.matrix((1, neles),
                                      extent=nonce + 'shock_normal_mag',
                                      tags={'align'})

        # shock_normal stores the *physical-space* unit vector for direct
        # visualisation in ParaView; the cascade converts the building block's
        # reference-space output via smats^T before writing to this matrix.
        eles.export_fields.append(ExportableField(
            name='shock-normal', shape=(eles.ndims,),
            getter=lambda: shock_normal.get().T
        ))
        eles.export_fields.append(ExportableField(
            name='shock-normal-mag', shape=(),
            getter=lambda: shock_normal_mag.get()[0]
        ))

        # Setup nodal/modal operator matrices
        invvdm, vdm_ef = self._build_operators(eles)

        if eles.basis.fpts_in_upts:
            m0 = None
        else:
            m0 = be.const_matrix(eles.basis.m0)

        # Jacobian data at upts -- "perfect information" for cascade development.
        # smats[i][j] = (J^-1 * det J)[i][j]; rcpdjac = 1/det J.  J^-1 itself is
        # `smats * rcpdjac`; J^-T is `smats^T * rcpdjac`; physical volume integ-
        # ration uses `det J = 1/rcpdjac`.  Per-element layout passed to the
        # kernel is [nupts][ndims*ndims] for smats and [nupts] for rcpdjac;
        # entry [q][i*ndims + j] is smats[i, j] at upt q.
        ndims = eles.ndims
        nupts = eles.nupts
        neles = eles.neles
        smats_np = eles.smat_at_np('upts').transpose(1, 0, 2, 3)
        smats_np = smats_np.reshape(nupts*ndims*ndims, neles)
        smats_upts = be.const_matrix(smats_np, tags={'align'})
        rcpdjac_upts = eles.rcpdjac_at('upts')

        # Reference-space gradient operator at upts (broadcast across elements).
        # Layout (ndims*nupts, nupts): row d*nupts + q, col i gives
        # d/d(xi_d) of nodal basis function i evaluated at upt q.  Applying to
        # a per-cell solution vector u_upts produces a stacked vector
        # [du/dxi_0 at all upts, du/dxi_1 at all upts, ...].
        #
        # Optionally fold in a modal projection that zeros out modes with
        # degree > shock-dir-p before differentiating.  The high modes carry
        # most of the polynomial Gibbs oscillation at a shock, so dropping
        # them gives a steadier direction estimate.  Default is the full
        # polynomial order (no projection).  Implemented as a precomputed
        # operator m4 @ P where P is the nodal->modal->truncate->nodal map,
        # so runtime cost is unchanged.
        shock_dir_p = cfg.getint('solver-entropy-filter', 'shock-dir-p',
                                 eles.basis.order)
        ub = eles.basis.ubasis
        sigma = np.array([1.0 if max(dd) <= shock_dir_p else 0.0
                          for dd in ub.degrees])
        proj = ub.vdm.T @ np.diag(sigma) @ ub.invvdm.T
        grad_op_np = eles.basis.m4 @ proj
        grad_op = be.const_matrix(grad_op_np, tags={'align'})

        # Build template arguments
        eftplargs = self._build_tplargs(eles, cfg, nfaces)

        # Register kernel factories on elements
        def local_entropy_kern(uin):
            return be.kernel(
                'entropylocal', tplargs=eftplargs, dims=[eles.neles],
                u=eles.scal_upts[uin], entmin_int=entmin, m0=m0
            )

        def entropy_filter_kern(uin):
            return be.kernel(
                'entropyfilter', tplargs=eftplargs, dims=[eles.neles],
                u=eles.scal_upts[uin], entmin_int=entmin, ef_filter=ef_filter,
                vdm=vdm_ef, invvdm=invvdm, m0=m0,
                mean_wts=eles.mean_wts,
                smats_upts=smats_upts, rcpdjac_upts=rcpdjac_upts,
                grad_op=grad_op,
                shock_normal=shock_normal,
                shock_normal_mag=shock_normal_mag
            )

        eles.kernels['local_entropy'] = local_entropy_kern
        eles.kernels['entropy_filter'] = entropy_filter_kern

    def _build_operators(self, eles):
        be = self._be

        invvdm = be.const_matrix(eles.basis.ubasis.invvdm.T)
        vdm = eles.basis.ubasis.vdm.T

        if not eles.basis.fpts_in_upts:
            vdmf = eles.basis.ubasis.vdm_at(eles.basis.fpts).T
            vdm = np.vstack([vdm, vdmf])

        return invvdm, be.const_matrix(vdm)

    def _build_tplargs(self, eles, cfg, nfaces):
        fpts_in_upts = eles.basis.fpts_in_upts
        nefpts = eles.nupts if fpts_in_upts else eles.nupts + eles.nfpts
        ub = eles.basis.ubasis

        # Reference-space face geometry for shock-normal detectors.
        face_ref_normals = [tuple(map(float, fn))
                            for _, _, fn in eles.basis.faces]
        face_ref_lengths = [_FACE_REF_MEASURE[k]
                            for k, _, _ in eles.basis.faces]
        cell_ref_volume = _CELL_REF_VOLUME[eles.basis.name]

        # Reference-space quadrature weights at upts (per-etype constants,
        # baked at template time).  Used by volume-integrated detectors.
        soln_pts_rule = cfg.get(f'solver-elements-{eles.basis.name}',
                                'soln-pts')
        upts_wts = list(map(float,
                            get_quadrule(eles.basis.name, soln_pts_rule,
                                         eles.nupts).wts))

        return {
            'ndims': eles.ndims, 'nupts': eles.nupts,
            'nfpts': eles.nfpts, 'nefpts': nefpts,
            'nvars': eles.nvars, 'nfaces': nfaces,
            'c': cfg.items_as('constants', float),
            'order': eles.basis.order,
            'fpts_in_upts': fpts_in_upts,
            'd_min': cfg.getfloat('solver-entropy-filter', 'd-min', 1e-6),
            'p_min': cfg.getfloat('solver-entropy-filter', 'p-min', 1e-6),
            'e_tol': cfg.getfloat('solver-entropy-filter', 'e-tol', 1e-6),
            'f_tol': cfg.getfloat('solver-entropy-filter', 'f-tol', 1e-4),
            'niters': cfg.getfloat('solver-entropy-filter', 'niters', 2),
            'cascade': cfg.get('solver-entropy-filter', 'cascade',
                               'legacy_nonlinear'),
            'ubdegs': [int(max(dd)) for dd in ub.degrees],
            # Reference-space face geometry for shock-normal detectors
            'face_ref_normals': face_ref_normals,
            'face_ref_lengths': face_ref_lengths,
            'cell_ref_volume': cell_ref_volume,
            'shock_normal_eps': 1e-14,
            # Volume-grad detector
            'upts_wts': upts_wts,
            'vg_weight_power': cfg.getint('solver-entropy-filter',
                                          'vg-weight-power', 1),
            'vg_field': cfg.get('solver-entropy-filter', 'vg-field',
                                'density'),
            # Face-grad detector (local divergence theorem with surface
            # quadrature over a chosen scalar field at all fpts).
            'facefpts': [list(map(int, ff)) for ff in eles.basis.facefpts],
            'fpts_wts': [float(w) for w in eles.basis.fpts_wts],
            'fg_field': cfg.get('solver-entropy-filter', 'fg-field',
                                'density'),
            # Directional convex filter
            'mode_ij': [tuple(int(d) for d in dd) for dd in ub.degrees],
            'dir_eps': cfg.getfloat('solver-entropy-filter', 'dir-eps',
                                    1e-3),
            # Which shock-normal detector to compile in.  See the building
            # blocks under entfilter/shock_normals/.
            'shock_normal_detector': cfg.get('solver-entropy-filter',
                                             'shock-normal',
                                             'volume_grad_density'),
        }

    def _setup_interfaces(self, system, int_inters, mpi_inters, bc_inters):
        be = self._be

        # Create interface views for entmin using system API
        iint_v, mpi_v, bc_v = system.make_field_views(
            self._entmin, layout='face', bc_layout='face-expand'
        )

        # Register MPI exchange for entropy face values
        system.register_mpi_exchange('ent_fpts', mpi_v)

        # Register comm_entropy kernels on internal/MPI interfaces
        def cent_kern(kname, intf, lhs, rhs):
            return lambda: be.kernel(kname, tplargs={}, dims=[intf.ninters],
                                     entmin_lhs=lhs, entmin_rhs=rhs)

        kprefix = 'pyfr.solvers.baseadvec.kernels'
        be.pointwise.register(f'{kprefix}.intcent')
        for i, (lhs, rhs) in zip(int_inters, iint_v):
            i.kernels['comm_entropy'] = cent_kern('intcent', i, lhs, rhs)

        be.pointwise.register(f'{kprefix}.mpicent')
        for m, (lhs, rhs) in zip(mpi_inters, mpi_v):
            m.kernels['comm_entropy'] = cent_kern('mpicent', m, lhs, rhs)

        # Register comm_entropy kernels on BC interfaces
        for b, lhs in zip(bc_inters, bc_v):
            b.kernels['comm_entropy'] = b.comm_entropy_kernel(lhs)

    def add_to_graph_pre_recv(self, g, k, m):
        # Post entropy MPI receives
        g.add_mpi_reqs(m['ent_fpts_recv'])

        # Run the entropy filter
        g.add_all(k['eles/entropy_filter'])

        # Pack and send entropy face values to neighbours
        g.add_all(k['mpiint/ent_fpts_pack'], deps=k['eles/entropy_filter'])
        for send, pack in zip(m['ent_fpts_send'], k['mpiint/ent_fpts_pack']):
            g.add_mpi_req(send, deps=[pack])

        # Compute common entropy at internal interfaces
        g.add_all(k['iint/comm_entropy'],
                  deps=k['eles/entropy_filter'] + k['mpiint/ent_fpts_pack'])

        # BC entropy needs the solution at flux points (disu)
        g.add_all(k['bcint/comm_entropy'], deps=k['eles/disu'])

    def add_to_graph_post_recv(self, g, k, deps):
        # Unpack MPI entropy data (may be empty when unpack is a no-op)
        g.add_all(k['mpiint/ent_fpts_unpack'])

        # Compute common entropy at MPI interfaces
        for c in k['mpiint/comm_entropy']:
            g.add(c, deps=deps(c, 'mpiint/ent_fpts_unpack'))

    def preproc_graphs(self, be, k, m, deps):
        # Graph: filter, compute local entropy, exchange, internal/BC entropy
        g_filter = be.graph()
        g_filter.add_mpi_reqs(m['ent_fpts_recv'])

        # Run entropy filter then compute local entropy minima
        g_filter.add_all(k['eles/entropy_filter'])
        g_filter.add_all(k['eles/local_entropy'],
                         deps=k['eles/entropy_filter'])

        # Interpolate to flux points (needed by bcint/comm_entropy)
        g_filter.add_all(k['eles/disu'], deps=k['eles/entropy_filter'])

        # Pack and send entropy values to neighbours
        g_filter.add_all(k['mpiint/ent_fpts_pack'],
                         deps=k['eles/local_entropy'])
        for send, pack in zip(m['ent_fpts_send'], k['mpiint/ent_fpts_pack']):
            g_filter.add_mpi_req(send, deps=[pack])

        # Compute common entropy minima at internal/boundary interfaces
        g_filter.add_all(k['iint/comm_entropy'], deps=k['eles/local_entropy'])
        g_filter.add_all(k['bcint/comm_entropy'],
                         deps=k['eles/local_entropy'] + k['eles/disu'])
        g_filter.commit()

        # Graph: MPI comm_entropy (if we have MPI interfaces)
        if k['mpiint/comm_entropy']:
            g_mpi_ent = be.graph()
            g_mpi_ent.add_all(k['mpiint/ent_fpts_unpack'])
            for c in k['mpiint/comm_entropy']:
                g_mpi_ent.add(c, deps=deps(c, 'mpiint/ent_fpts_unpack'))

            g_mpi_ent.commit()
            return g_filter, g_mpi_ent

        return g_filter,

    def postproc(self, be, k):
        be.run_kernels(k['eles/entropy_filter'])
