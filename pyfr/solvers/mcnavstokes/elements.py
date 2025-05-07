import numpy as np

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.mceuler.elements import BaseMCFluidElements
from pyfr.multicomp.mcfluid import MCFluid


class MCNavierStokesElements(BaseMCFluidElements,
                             BaseAdvectionDiffusionElements):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.mcfluid = MCFluid(self.cfg, justTherm=False)

    @staticmethod
    def grad_con_to_pri(cons, grad_cons, cfg):
        fluid = MCFluid(cfg, justTherm=True)
        return fluid.diff_con_to_pri(cons, grad_cons)

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        kprefix = 'pyfr.solvers.mcnavstokes.kernels'
        self._be.pointwise.register(f'{kprefix}.tflux')

        # Handle shock capturing
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        # Viscous Sponge
        visc_sponge = 'viscous-sponge' in self.cfg.sections()

        # Template parameters for the flux kernels
        consts = self.cfg.items_as('constants', float)
        consts |= self.mcfluid.consts
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'c': consts,
            'eos': self.mcfluid.eos,
            'trans': self.mcfluid.trans,
            'mixing_rule': self.mcfluid.mixing_rule,
            'jac_exprs': self.basis.jac_exprs,
            'shock_capturing': shock_capturing,
            'visc_sponge': visc_sponge,
        }

        if visc_sponge:
            start_point = np.asarray([float(i) for i in self.cfg.get('viscous-sponge','start').split(',')])
            end_point = np.asarray([float(i) for i in self.cfg.get('viscous-sponge','end').split(',')])
            mult = self.cfg.getfloat('viscous-sponge','mult')
            profile = self.cfg.get('viscous-sponge','profile', default='linear')
            ploc = self.ploc_at_np('upts')
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

        # Helpers
        tdisf = []
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat
        av = self.artvisc

        # Gradient + flux kernel fusion
        if self.grad_fusion:
            if c in r:
                tdisf.append(lambda uin: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'curved-fused'},
                    dims=[self.nupts, r[c]], u=s(self.scal_upts[uin], c),
                    artvisc=s(av, c), f=s(self._vect_upts, c),
                    gradu=s(self._grad_upts, c),
                    rcpdjac=self.rcpdjac_at('upts', 'curved'),
                    smats=self.curved_smat_at('upts'),
                    extrns=self._external_args,
                    **self._external_vals
                ))
            if l in r:
                tdisf.append(lambda uin: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'linear-fused'},
                    dims=[self.nupts, r[l]], u=s(self.scal_upts[uin], l),
                    artvisc=s(av, l), f=s(self._vect_upts, l),
                    gradu=s(self._grad_upts, l),
                    verts=self.ploc_at('linspts', l), upts=self.upts,
                    extrns=self._external_args,
                    **self._external_vals
                ))

            def tdisf_k(uin):
                return self._make_sliced_kernel(k(uin) for k in tdisf)

            self.kernels['tdisf_fused'] = tdisf_k
        # No gradient + flux kernel fusion, with flux-AA
        elif 'flux' in self.antialias:
            if c in r:
                tdisf.append(lambda: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'curved'},
                    dims=[self.nqpts, r[c]], u=s(self._scal_qpts, c),
                    f=s(self._vect_qpts, c), artvisc=s(av, c),
                    smats=self.curved_smat_at('qpts'),
                    extrns=self._external_args,
                    **self._external_vals
                ))
            if l in r:
                tdisf.append(lambda: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'linear'},
                    dims=[self.nqpts, r[l]], u=s(self._scal_qpts, l),
                    f=s(self._vect_qpts, l), artvisc=s(av, l),
                    verts=self.ploc_at('linspts', l), upts=self.qpts,
                    extrns=self._external_args,
                    **self._external_vals
                ))

            def tdisf_k():
                return self._make_sliced_kernel(k() for k in tdisf)

            self.kernels['tdisf'] = tdisf_k
        # No gradient + flux kernel fusion, no flux-AA
        else:
            if c in r:
                tdisf.append(lambda uin: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'curved'},
                    dims=[self.nupts, r[c]], u=s(self.scal_upts[uin], c),
                    f=s(self._vect_upts, c), artvisc=s(av, c),
                    smats=self.curved_smat_at('upts'),
                    extrns=self._external_args,
                    **self._external_vals
                ))
            if l in r:
                tdisf.append(lambda uin: self._be.kernel(
                    'tflux', tplargs=tplargs | {'ktype': 'linear'},
                    dims=[self.nupts, r[l]], u=s(self.scal_upts[uin], l),
                    f=s(self._vect_upts, l), artvisc=s(av, l),
                    verts=self.ploc_at('linspts', l), upts=self.upts,
                    extrns=self._external_args,
                    **self._external_vals
                ))

            def tdisf_k(uin):
                return self._make_sliced_kernel(k(uin) for k in tdisf)

            self.kernels['tdisf'] = tdisf_k