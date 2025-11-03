import numpy as np

from pyfr.util import first


def _get_inter_objs(interside, getter, elemap):
    # Map from element type to view mat getter
    emap = {type: getattr(ele, getter) for type, ele in elemap.items()}

    # Get the data from the interface
    return [emap[type](eidx, fidx) for type, eidx, fidx in interside]


class BaseInters:
    def __init__(self, be, lhs, elemap, cfg):
        self._be = be
        self.elemap = elemap
        self.cfg = cfg

        # Get the number of dimensions and variables
        self.ndims = first(elemap.values()).ndims
        self.nvars = first(elemap.values()).nvars

        # Get the number of interfaces
        self.ninters = len(lhs)

        # Compute the total number of interface flux points
        self.ninterfpts = sum(elemap[etype].nfacefpts[fidx]
                              for etype, eidx, fidx in lhs)

        # By default do not permute any of the interface arrays
        self._perm = Ellipsis

        # Kernel constants
        self.c = cfg.items_as('constants', float)

        # Kernels and MPI requests we provide
        self.kernels = {}
        self.mpireqs = {}

        # Global kernel arguments
        self._external_args = {}
        self._external_vals = {}

        # Viscous Sponge
        self.visc_sponge = 'viscous-sponge' in self.cfg.sections()

        if self.visc_sponge:
            import numpy as np
            start_point = np.asarray([float(i) for i in self.cfg.get('viscous-sponge','start').split(',')])
            end_point = np.asarray([float(i) for i in self.cfg.get('viscous-sponge','end').split(',')])
            mult = self.cfg.getfloat('viscous-sponge','mult')
            profile = self.cfg.get('viscous-sponge','profile', default='linear')
            ploc = self._const_mat(lhs, 'get_ploc_for_inter').get()
            # Calculate the direction vector from start to end
            direction = end_point - start_point
            direction_norm = np.linalg.norm(direction)

            if direction_norm < 1e-10:  # Start and end points are too close
                raise ValueError("Sponge start and end points are too close or identical")

            # Normalize the direction vector
            direction_unit = direction / direction_norm

            # For each point, calculate its projection onto the line from start to end
            # First, vector from start to each point
            vectors_from_start = ploc - start_point[np.newaxis,:,np.newaxis]

            # Project these vectors onto the direction unit vector
            projections = np.einsum('ijk,j->ik', vectors_from_start, direction_unit)

            # Normalize projections to get a parameter t between 0 and 1
            # where t=0 at start_point and t=1 at end_point
            t = projections / direction_norm

            # Apply the growth function
            if profile.lower() == 'linear':
                # Linear growth from 1 to mag
                scaled = np.clip(t, 0, 1)
            elif profile.lower() == 'quadratic':
                # Quadratic growth from 1 to mag
                scaled = np.clip(t, 0, 1)**2
            elif profile.lower() == 'tanh':
                # Scale t from [0,1] to [-1.5,1.5]
                # Apply tanh and then rescale to ensure we go exactly from 1 to mag
                scaled = (np.tanh(3 * (t - 0.5)) - np.tanh(-1.5)) / (np.tanh(1.5) - np.tanh(-1.5))
            else:
                raise ValueError("Sponge growth type must be 'linear', 'quadratic', or 'tanh'")
            multipliers = 1.0 + (mult - 1.0) * scaled

            # Set multiplier to 1 for all points before start point (t < 0)
            # and to mag for all points after end point (t > 1)
            multipliers = np.where(t < 0, 1.0, multipliers)
            multipliers = np.where(t > 1, mult, multipliers)

            self._set_external('sponge_mult',
                               'in fpdtype_t',
                               self._be.const_matrix(multipliers))

    def _set_external(self, name, spec, value=None):
        self._external_args[name] = spec

        if value is not None:
            self._external_vals[name] = value

    def _const_mat(self, inter, meth):
        m = _get_inter_objs(inter, meth, self.elemap)

        # Swizzle the dimensions and permute
        m = np.concatenate(m)
        m = np.atleast_2d(m.T)
        m = m[:, self._perm]

        return self._be.const_matrix(m)

    def _ewise_const_mat(self, inter, meth):
        m = _get_inter_objs(inter, meth, self.elemap)

        # Swizzle the dimensions
        m = np.array(m)
        m = np.moveaxis(m, 0, -1)

        return self._be.const_matrix(m)

    def _get_perm_for_view(self, inter, meth):
        vm = _get_inter_objs(inter, meth, self.elemap)
        vm = [np.concatenate(m) for m in zip(*vm)]
        mm = self._be.view(*vm, vshape=()).mapping.get()

        return np.argsort(mm[0])

    def _view(self, inter, meth, vshape=(), with_perm=True):
        vm = _get_inter_objs(inter, meth, self.elemap)
        perm = self._perm if with_perm else Ellipsis
        vm = [np.concatenate(m)[perm] for m in zip(*vm)]
        return self._be.view(*vm, vshape=vshape)

    def _scal_view(self, inter, meth):
        return self._view(inter, meth, (self.nvars,))

    def _vect_view(self, inter, meth):
        return self._view(inter, meth, (self.ndims, self.nvars))

    def _scal_upts_view(self, inter, meth):
        nupts = first(self.elemap.values()).basis.nupts
        return self._view(inter, meth, (nupts, self.nvars), with_perm=False)

    def _scal_fpts_view(self, inter, meth):
        basis = first(self.elemap.values()).basis
        vshape = (basis.nfpts, self.nvars)
        with_perm = False
        return self._view(inter, meth, vshape=vshape, with_perm=with_perm)

    def _grad_upts_view(self, inter, meth):
        basis = first(self.elemap.values()).basis
        vshape = (self.ndims*basis.nupts, self.nvars)
        with_perm = False
        return self._view(inter, meth, vshape=vshape, with_perm=with_perm)

    def _vect_fpts_view(self, inter, meth):
        basis = first(self.elemap.values()).basis
        vshape = (self.ndims*basis.nfpts, self.nvars)
        with_perm = False
        return self._view(inter, meth, vshape=vshape, with_perm=with_perm)

    def _xchg_view(self, inter, meth, vshape=(), with_perm=True):
        vm = _get_inter_objs(inter, meth, self.elemap)
        perm = self._perm if with_perm else Ellipsis
        vm = [np.concatenate(m)[perm] for m in zip(*vm)]
        return self._be.xchg_view(*vm, vshape=vshape)

    def _scal_xchg_view(self, inter, meth):
        return self._xchg_view(inter, meth, (self.nvars,))

    def _vect_xchg_view(self, inter, meth):
        return self._xchg_view(inter, meth, (self.ndims, self.nvars))
    
    def setup(self, sdata):
        pass

    @classmethod
    def serialisefn(cls, iface, prefix, srl):
        pass
