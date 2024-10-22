import numpy as np

from pyfr.mpiutil import get_comm_rank_root
from pyfr.points import PointSampler
from pyfr.polys import TriPolyBasis
from pyfr.shapes import TriShape
from pyfr.writers.vtk.base import BaseVTKWriter, interpolate_pts


class VTKSTLWriter(BaseVTKWriter):
    type = 'stl'

    def __init__(self, meshf, solnf, stlrgns, **kwargs):
        # Disable high-order output and (by default) subdivision
        kwargs['order'] = None
        kwargs.setdefault('divisor', 1)

        super().__init__(meshf, solnf, **kwargs)

        if self.ndims != 3:
            raise RuntimeError('STL export only supported for 3D grids')

        mesh, soln = self.mesh, self.soln
        comm, rank, root = get_comm_rank_root()

        # Read and merge the STL surfaces from the mesh
        pts = np.vstack([mesh.raw[f'regions/stl/{s}'][:, 1:] for s in stlrgns])

        # Subdivide the mesh
        pts = self._subdivide_pts(pts.swapaxes(0, 1), self.divisor)

        # Determine the unique vertices
        ppts, pinv = np.unique(pts.reshape(-1, 3), axis=0, return_inverse=True)

        # Create a configure a point sampler
        sampler = PointSampler(mesh, ppts)
        sampler.configure_with_cfg_nvars(soln['config'],
                                         len(self._soln_fields))

        # Perform the sampling
        samps = sampler.sample([soln[etype] for etype in mesh.eidxs])

        # If we are the root rank then write out the triangle list
        if rank == root:
            # Pre-process the samples
            samps = samps.swapaxes(0, 1).astype(self.dtype)
            samps = self._pre_proc_fields(samps)

            # Unpack and reshape
            soln = samps[:, pinv].reshape(-1, *pts.shape[:2]).swapaxes(0, 1)

            self.einfo = [('tri', pts.shape[1])]
            self._stl_info = pts, soln
        else:
            self.einfo = []

    def _prepare_pts(self, etype):
        return *self._stl_info, None, None

    def _subdivide_pts(self, pts, order):
        basis = TriPolyBasis(2, TriShape.std_ele(1))
        op = basis.nodal_basis_at(TriShape.std_ele(order))

        return interpolate_pts(op, pts)
