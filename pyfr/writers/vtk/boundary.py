from collections import defaultdict

import numpy as np

from pyfr.shapes import BaseShape
from pyfr.util import subclass_where
from pyfr.writers.vtk.base import BaseVTKWriter, interpolate_pts


def _search(a, v):
    idx = np.argsort(a)
    return idx[np.searchsorted(a, v, sorter=idx)]


class VTKBoundaryWriter(BaseVTKWriter):
    type = 'boundary'
    output_curved = True
    output_partition = True

    def __init__(self, meshf, solnf, boundaries, **kwargs):
        super().__init__(meshf, solnf, **kwargs)

        if self.ndims != 3:
            raise RuntimeError('Boundary export only supported for 3D grids')

        ecount = defaultdict(int)
        self._surface_info = defaultdict(list)

        rmesh, smesh = self.reader.mesh, self.mesh
        cidxs = [smesh.codec.index(f'bc/{b}') for b in boundaries]

        for etype, einfo in self.reader.eles.items():
            # See which of our faces are on the selected boundaries
            mask = np.isin(einfo['faces']['cidx'], cidxs)
            eoff, fidx = mask.nonzero()

            # Handle the case where the solution is subset
            if smesh.subset:
                eidxs = rmesh.eidxs[etype]
                beidx = eidxs[mask.any(axis=1)]
                if not np.isin(beidx, smesh.eidxs[etype]).all():
                    raise ValueError('Output boundaries not present in subset '
                                     'solution')

                eoff = _search(smesh.eidxs[etype], eidxs[eoff])

            # Obtain the associated surface info
            for stype, *info in self._get_surface_info(etype, eoff, fidx):
                ecount[stype] += len(info[-1])
                self._surface_info[stype].append((etype, *info))

        self.einfo = list(ecount.items())

    def _get_surface_info(self, etype, eoff, fidx):
        info, idxs = {}, defaultdict(list)

        # Get the shape associated with our element type
        shapecls = subclass_where(BaseShape, name=etype)
        shape = shapecls(len(self.mesh.spts[etype]), self.cfg)

        for e, f in zip(eoff, fidx):
            if f not in info:
                itype, proj, norm = shape.faces[f]
                ishapecls = subclass_where(BaseShape, name=itype)

                # Obtain the visualisation points on this face
                svpts = np.array(ishapecls.std_ele(self.etypes_div[itype]))
                svpts = np.vstack(np.broadcast_arrays(*proj(*svpts.T))).T

                if self.ho_output:
                    svpts = svpts[self._nodemaps[itype, len(svpts)]]

                mesh_op = shape.sbasis.nodal_basis_at(svpts)
                soln_op = shape.ubasis.nodal_basis_at(svpts)

                info[f] = (itype, mesh_op, soln_op)

            idxs[f].append(e)

        return [(*info[f], idxs[f]) for f in info]

    def _prepare_pts(self, itype):
        vspts, vsoln, curved, part = [], [], [], []

        for etype, mesh_op, soln_op, idxs in self._surface_info[itype]:
            spts = self.mesh.spts[etype][:, idxs]
            soln = self.soln[etype][..., idxs]
            soln = soln.swapaxes(0, 1).astype(self.dtype)

            # Pre-process the solution
            soln = self._pre_proc_fields(soln).swapaxes(0, 1)

            vspts.append(interpolate_pts(mesh_op, spts))
            vsoln.append(interpolate_pts(soln_op, soln))
            curved.append(self.mesh.spts_curved[etype][idxs])
            part.append(self.soln[f'{etype}-parts'][idxs])

        return (np.hstack(vspts), np.dstack(vsoln),
                np.hstack(curved), np.hstack(part))