from math import prod

from pyfr.backends.base.generator import BaseGPUKernelGenerator
from pyfr.backends.base.ikpgenerator import GPUIKPKernelGeneratorMixin


class CUDAKernelGenerator(BaseGPUKernelGenerator, GPUIKPKernelGeneratorMixin):
    _lid = ('threadIdx.x', 'threadIdx.y')
    _gid = 'ixdtype_t(blockIdx.x)*blockDim.x + threadIdx.x'
    _shared_prfx = '__shared__'
    _shared_sync = '__syncthreads()'

    def _render_spec(self):
        res = '__restrict__'

        # We first need the argument list; starting with the dimensions
        kargs = [f'ixdtype_t {d}' for d in self._dims]

        # Now add any scalar arguments
        kargs.extend(f'{sa.dtype} {sa.name}' for sa in self.scalargs)

        # Finally, add the vector arguments
        for va in self.vectargs:
            if va.intent == 'in':
                kargs.append(f'const {va.dtype}* {res} {va.name}_v')
            else:
                kargs.append(f'{va.dtype}* {res} {va.name}_v')

            # Views
            if va.isview:
                kargs.append(f'const ixdtype_t* {res} {va.name}_vix')

                if va.ncdim == 2:
                    kargs.append(f'const ixdtype_t* {res} {va.name}_vrstri')
            # Arrays
            elif self.needs_ldim(va):
                kargs.append(f'ixdtype_t ld{va.name}')

        # Determine the launch bounds for the kernel
        nthrds = prod(self.block1d if self.ndim == 1 else self.block2d)
        kattrs = f'__global__ __launch_bounds__({nthrds})'

        return '{0} void {1}({2})'.format(kattrs, self.name, ', '.join(kargs))

    def _render_body_preamble(self, body):
        """
        Transform kernel body and generate preamble.

        Calls base class for standard dereferencing, then performs GPU IKP
        transformations if enabled.
        """
        # Standard dereferencing (base class)
        body, preamble = super()._render_body_preamble(body)

        # Apply GPU IKP transformations if enabled
        if self.ikp:
            body, preamble = self._ikp_render_body_preamble(body, preamble)

        return body, preamble

    def render(self):
        """
        Render the complete kernel.

        For IKP mode, uses cooperative thread groups where multiple threads
        work on the same element.
        """
        spec = self._render_spec()

        # Use IKP-specific global ID calculation if IKP is enabled
        gid = self._ikp_gid() if self.ikp else self._gid

        return f'''{spec}
            {{
                ixdtype_t _x = {gid};
                #define X_IDX (_x)
                #define X_IDX_AOSOA(v, nv) SOA_IX(X_IDX, v, nv)
                #define BCAST_BLK(r, c, ld)  c
                {self.preamble}
                {{
                    {self.body}
                }}
                #undef X_IDX
                #undef X_IDX_AOSOA
                #undef BCAST_BLK
            }}'''
