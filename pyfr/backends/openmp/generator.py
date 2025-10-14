from math import prod
import re

from pyfr.backends.base.generator import BaseKernelGenerator
from pyfr.backends.base.ikpgenerator import IKPKernelGeneratorMixin


class OpenMPKernelGenerator(IKPKernelGeneratorMixin, BaseKernelGenerator):
    """
    OpenMP kernel generator with optional IKP (Inner-Kernel Parallelism) support.

    Standard mode: SIMD parallelism across elements (SOA_SZ chunks)
    IKP mode: Sequential element loop with cache-blocking

    IKP Strategy:
    - Process BLK_SZ elements sequentially to maximize cache reuse
    - Hoist arrays to block-level: arr[N] -> arr[BLK_SZ][N]
    - Reference with ELEM_IDX: arr[i] -> arr[ELEM_IDX][i] where ELEM_IDX = _elem
    - Future: Replace unrolled GEMV with batched libxsmm calls

    Note: IKP is CPU-specific. GPU backends will use thread cooperation instead.
    """

    # ========== IKP Mixin Implementation (OpenMP-Specific) ==========

    def _ikp_transform_const_decl(self, dtype, name, dims, initializer):
        """
        Transform constant array declaration for OpenMP IKP.

        Constants are shared across all elements (no transformation needed).
        """
        return f'    const {dtype} {name}{dims} = {initializer};'

    def _ikp_transform_local_decl(self, dtype, name, size):
        """
        Transform local array declaration for OpenMP IKP.

        arr[N] -> arr[N][BLK_SZ] (column-major for libxsmm compatibility)
        """
        return f'    {dtype} {name}[{size}][BLK_SZ];'

    def _ikp_transform_array_ref(self, arr_name, body):
        """
        Transform local array references for OpenMP IKP.

        arr[i] -> arr[i][X_IDX] (column-major layout for libxsmm)

        This layout matches libxsmm's expected transposed format where
        BLK_SZ elements are stored contiguously for each array index.
        """
        arr_pattern = rf'\b{arr_name}\[([^\]]+)\]'
        return re.sub(arr_pattern, rf'{arr_name}[\1][X_IDX]', body)

    def _ikp_transform_kernel_args(self, body):
        """
        Transform kernel argument references for OpenMP IKP.

        No transformation needed - X_IDX and X_IDX_AOSOA already refer to
        the correct element index for parallel execution.
        """
        # No changes needed - keep using X_IDX and X_IDX_AOSOA as-is
        return body

    def _ikp_elem_idx_macro(self):
        """
        Return ELEM_IDX macro definition for OpenMP.

        Not needed anymore - we use X_IDX directly.
        """
        return ''

    def _ikp_wrap_body(self, body, nelem_expr):
        """
        Wrap IKP body in standard _xi/_xj loop structure for parallel execution.

        Splits body into prep/gemm/proc phases and structures them correctly:
        1. Prep (extract data) - INSIDE element loops
        2. GEMM (libxsmm) - OUTSIDE loops (batched over all elements)
        3. Proc (process results) - INSIDE element loops
        """
        # Remove IKP_LOOP markers
        body = re.sub(r'// IKP_LOOP_BEGIN\n', '', body)
        body = re.sub(r'// IKP_LOOP_END', '', body)

        # Split body into prep and proc sections at the prep_end boundary
        # Pattern: everything before the gemv block opening brace
        prep_match = re.search(r'(.*?)\s*\{\s*// PYFR_IKP_MARKER:', body, flags=re.DOTALL)
        prep_body = prep_match.group(1) if prep_match else body

        # Extract GEMM section (from IKP_MARKER to gemm_end)
        gemm_match = re.search(
            r'(// PYFR_IKP_MARKER:.*?// PYFR_IKP_PHASE_BOUNDARY: gemm_end\n)',
            body, flags=re.DOTALL
        )
        gemm_section = gemm_match.group(1) if gemm_match else ''

        # Proc section: everything after gemm_end, skip closing brace
        proc_match = re.search(
            r'// PYFR_IKP_PHASE_BOUNDARY: gemm_end\n\s*\}\s*(.*)',
            body, flags=re.DOTALL
        )
        proc_body = proc_match.group(1) if proc_match else ''

        # Create loop wrappers - prep and proc need different _xi handling
        def make_prep_loop(content):
            if nelem_expr == 'BLK_SZ':
                # Core path: full block
                return f'''
                for (int _xi = 0; _xi < BLK_SZ; _xi += SOA_SZ)
                {{
                    #pragma omp simd
                    for (int _xj = 0; _xj < SOA_SZ; _xj++)
                    {{
                        {content}
                    }}
                }}'''
            else:
                # Clean path: partial block - declare _xi here
                return f'''
                int _xi = 0;
                #pragma omp simd
                for (int _xj = 0; _xj < _nx % BLK_SZ; _xj++)
                {{
                    {content}
                }}'''

        def make_proc_loop(content):
            if nelem_expr == 'BLK_SZ':
                # Core path: full block
                return f'''
                for (int _xi = 0; _xi < BLK_SZ; _xi += SOA_SZ)
                {{
                    #pragma omp simd
                    for (int _xj = 0; _xj < SOA_SZ; _xj++)
                    {{
                        {content}
                    }}
                }}'''
            else:
                # Clean path: reuse _xi from prep loop
                return f'''
                #pragma omp simd
                for (int _xj = 0; _xj < _nx % BLK_SZ; _xj++)
                {{
                    {content}
                }}'''

        # Assemble: prep loop -> gemm -> proc loop
        # Note: GEMM goes between the two loops, NOT nested inside
        return make_prep_loop(prep_body) + '\n' + gemm_section + '\n' + make_proc_loop(proc_body)

    # ========== IKP Transformation Pipeline ==========

    def _render_body_preamble(self, body):
        """
        Override to hoist declarations for IKP mode.

        Calls base class for standard dereferencing, then hoists IKP arrays.
        """
        # Standard dereferencing (base class)
        body, preamble = super()._render_body_preamble(body)

        if self.ikp:
            # Find declarations using mixin methods
            const_arrays, local_arrays = self._ikp_find_declarations(body)

            # Track local array names for reference transformation
            self._ikp_local_arrays = [name for _, name, _ in local_arrays]

            # Transform declarations
            declarations = []

            # Transform constant arrays (shared across elements)
            for dtype, name, dims, initializer in const_arrays:
                decl = self._ikp_transform_const_decl(dtype, name, dims, initializer)
                declarations.append(decl)

            # Transform local arrays (per-element storage)
            for dtype, name, size in local_arrays:
                decl = self._ikp_transform_local_decl(dtype, name, size)
                declarations.append(decl)

            # Remove declarations from body (they're now in preamble)
            body = self._ikp_remove_declarations(body)

            # Add transformed declarations to preamble
            if declarations:
                preamble = '\n'.join(declarations) + '\n' + preamble

        return body, preamble

    def _transform_ikp_body(self, body):
        """
        Transform body for IKP: update array refs and kernel arg refs.

        Only transforms IKP_LOOP sections (preserves IKP_GEMM for future batching).
        """
        local_arrays = getattr(self, '_ikp_local_arrays', [])

        def transform_section(match):
            section_code = match.group(1)

            # Transform local array references
            for arr_name in local_arrays:
                section_code = self._ikp_transform_array_ref(arr_name, section_code)

            # Transform kernel argument references
            section_code = self._ikp_transform_kernel_args(section_code)

            return f'// IKP_LOOP_BEGIN\n{section_code}// IKP_LOOP_END'

        return re.sub(r'// IKP_LOOP_BEGIN\n(.*?)// IKP_LOOP_END',
                      transform_section, body, flags=re.DOTALL)

    # ========== Standard OpenMP Generator Methods ==========

    def render(self):
        kargdefn, kargassn = self._render_args('args')

        if self.ikp:
            # IKP mode: sequential element loop for cache blocking
            transformed_body = self._transform_ikp_body(self.body)
            core = self._ikp_wrap_body(transformed_body, 'BLK_SZ')
            clean = self._ikp_wrap_body(transformed_body, '(_nx % BLK_SZ)')
        elif self.ndim == 1:
            # Standard 1D: SIMD across elements
            core = f'''
                for (int _xi = 0; _xi < BLK_SZ; _xi += SOA_SZ)
                {{
                    #pragma omp simd
                    for (int _xj = 0; _xj < SOA_SZ; _xj++)
                    {{
                        {self.body}
                    }}
                }}'''
            clean = f'''
                int _xi = 0;
                #pragma omp simd
                for (int _xj = 0; _xj < _nx % BLK_SZ; _xj++)
                {{
                    {self.body}
                }}'''
        else:
            # Standard 2D: SIMD across elements
            core = f'''
                for (ixdtype_t _y = 0; _y < _ny; _y++)
                {{
                    for (int _xi = 0; _xi < BLK_SZ; _xi += SOA_SZ)
                    {{
                        #pragma omp simd
                        for (int _xj = 0; _xj < SOA_SZ; _xj++)
                        {{
                            {self.body}
                        }}
                    }}
                }}'''
            clean = f'''
                for (ixdtype_t _y = 0, _xi = 0; _y < _ny; _y++)
                {{
                    #pragma omp simd
                    for (int _xj = 0; _xj < _nx % BLK_SZ; _xj++)
                    {{
                        {self.body}
                    }}
                }}'''

        # Define macros (standard + IKP if needed)
        macros_def = '''#define X_IDX (_xi + _xj)
                #define X_IDX_AOSOA(v, nv)\
                    ((_xi/SOA_SZ*(nv) + (v))*SOA_SZ + _xj)
                #define BCAST_BLK(r, c, ld) ((c) % (ld) + ((c) / (ld))*(ld)*r)'''
        macros_undef = '''#undef X_IDX
                #undef X_IDX_AOSOA
                #undef BCAST_BLK'''

        # IKP mode uses X_IDX directly, no additional macros needed

        return f'''
            struct {self.name}_kargs {{ {kargdefn}; }};
            void {self.name}(ixdtype_t _ib,
                             const struct {self.name}_kargs *args,
                             int _disp_mask)
            {{
                {kargassn};
                {macros_def}
                {self.preamble}
                if (_nx - _ib*BLK_SZ >= BLK_SZ)
                {{
                    {core}
                }}
                else
                {{
                    {clean}
                }}
                {macros_undef}
            }}'''

    def ldim_size(self, name, factor=1):
        return f'{factor}*BLK_SZ' if factor > 1 else 'BLK_SZ'

    def needs_ldim(self, arg):
        return False

    def _displace_arg(self, arg):
        if arg.isview:
            return None
        elif self.ndim == 1:
            # Vector
            if arg.ncdim == 0 or arg.ismpi:
                return '_ib*BLK_SZ'
            # 2D broadcast vector
            elif arg.isbroadcast:
                return None
            # Stacked vector:
            else:
                return f'_ib*BLK_SZ*{prod(arg.cdims)}'
        else:
            # 2D broadcast vector or row broadcast matrix
            if arg.isbroadcast or arg.isbroadcastr:
                return None
            # Column broadcast matrix
            elif arg.isbroadcastc:
                return f'_ib*BLK_SZ*{prod(arg.cdims)}'
            # Matrix
            else:
                return f'_ib*BLK_SZ*{prod(arg.cdims)}*_ny'

    def _render_args(self, argn):
        # We first need the argument list; starting with the dimensions
        kargs = [('ixdtype_t', d, None, None) for d in self._dims]

        # Now add any scalar arguments
        kargs.extend((sa.dtype, sa.name, None, None) for sa in self.scalargs)

        # Finally, add the vector arguments
        for va in self.vectargs:
            da = self._displace_arg(va)
            mi = len(kargs) if da else None

            if va.intent == 'in':
                kargs.append((f'const {va.dtype}*', f'{va.name}_v', da, mi))
            else:
                kargs.append((f'{va.dtype}*', f'{va.name}_v', da, mi))

            # Views
            if va.isview:
                kargs.append(('const ixdtype_t*', f'{va.name}_vix',
                              '_ib*BLK_SZ', None))

                if va.ncdim == 2:
                    kargs.append(('const ixdtype_t*', f'{va.name}_vrstri',
                                  '_ib*BLK_SZ', None))

        # Argument definitions and assignments
        kargdefn, kargassn = [], []
        for dtype, name, disp, midx in kargs:
            assn = f'{dtype} {name} = {argn}->{name}'

            # Handle displacement and potential masking thereof
            if disp and midx is not None:
                assn += f' + ((_disp_mask & {1 << midx}) ? 0 : {disp})'
            elif disp:
                assn += f' + {disp}'

            kargdefn.append(f'{dtype} {name}')
            kargassn.append(assn)

        return ';\n'.join(kargdefn), ';\n'.join(kargassn)
