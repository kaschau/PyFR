import re

from pyfr.backends.base.generator import BaseGPUKernelGenerator
from pyfr.backends.base.ikpgenerator import GPUIKPKernelGeneratorMixin


class MetalKernelGenerator(BaseGPUKernelGenerator, GPUIKPKernelGeneratorMixin):
    _lid = ('_tpitg.x', '_tpitg.y')
    _gid = '_tpig.x'
    _shared_prfx = 'threadgroup'
    _shared_sync = 'threadgroup_barrier(mem_flags::mem_threadgroup)'

    def _render_spec(self):
        # We first need the argument list; starting with the dimensions
        kargs = [f'constant ixdtype_t& {d}' for d in self._dims]

        # Now add any scalar arguments
        kargs.extend(f'constant {sa.dtype}& {sa.name}' for sa in self.scalargs)

        # Then vector arguments
        for va in self.vectargs:
            if va.intent == 'in':
                kargs.append(f'device const {va.dtype}* {va.name}_v')
            else:
                kargs.append(f'device {va.dtype}* {va.name}_v')

            # Views
            if va.isview:
                kargs.append(f'device const ixdtype_t* {va.name}_vix')

                if va.ncdim == 2:
                    kargs.append(f'device const ixdtype_t* {va.name}_vrstri')
            # Arrays
            elif self.needs_ldim(va):
                kargs.append(f'device const ixdtype_t& ld{va.name}')

        # Finally, the attribute arguments
        kargs.append('uint2 _tpig [[thread_position_in_grid]]')
        # IKP kernels always need _tpitg for thread cooperation
        if self.ikp or re.search(r'\b_tpitg\b', self.preamble):
            kargs.append('uint2 _tpitg [[thread_position_in_threadgroup]]')

        return 'kernel void {0}({1})'.format(self.name, ', '.join(kargs))

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
