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

        arr[N] -> arr[N*BLK_SZ] (flat array like split version)
        """
        return f'    {dtype} {name}[BLK_SZ*{size}];'

    def _ikp_transform_array_ref(self, var_name, body):
        """
        Transform local variable references for OpenMP IKP.

        Arrays: arr[i] -> arr[i*BLK_SZ + X_IDX] (flat indexing like split)
        Scalars: scalar -> scalar[X_IDX] (elevated to array[BLK_SZ])

        Using flat indexing matches the split version and helps the compiler
        recognize the access pattern for better vectorization.
        """
        # First try array pattern (has brackets)
        arr_pattern = rf'\b{var_name}\[([^\]]+)\]'
        if re.search(arr_pattern, body):
            # It's an array reference - transform arr[i] -> arr[i*BLK_SZ + X_IDX]
            return re.sub(arr_pattern, rf'{var_name}[(\1)*BLK_SZ + X_IDX]', body)
        else:
            # It's a scalar reference - transform scalar -> scalar[X_IDX]
            # But be careful not to transform the declaration itself
            scalar_pattern = rf'\b{var_name}\b(?!\s*\[)'
            return re.sub(scalar_pattern, rf'{var_name}[X_IDX]', body)

    def _ikp_wrap_body(self, body, nelem_expr):
        """
        Wrap IKP body in standard _xi/_xj loop structure for parallel execution.

        Splits body into alternating prep/interruption sections:
        - Prep sections (per-element code) - INSIDE element loops
        - Interruption sections (batched ops) - OUTSIDE loops
        Handles N interruptions automatically.
        """
        # Remove IKP_LOOP markers
        body = re.sub(r'// IKP_LOOP_BEGIN\n', '', body)
        body = re.sub(r'// IKP_LOOP_END', '', body)

        # Split body on interruption markers
        # Pattern: PYFR_IKP_INTERRUPTION_START ... PYFR_IKP_INTERRUPTION_END
        sections = []
        remaining = body

        while True:
            # Find next interruption
            match = re.search(
                r'(.*?)// PYFR_IKP_INTERRUPTION_START\n(.*?)// PYFR_IKP_INTERRUPTION_END',
                remaining, flags=re.DOTALL
            )

            if not match:
                # No more interruptions, rest is final proc section
                if remaining.strip():
                    sections.append(('prep', remaining))
                break

            # Extract prep section before interruption
            prep_section = match.group(1)
            if prep_section.strip():
                sections.append(('prep', prep_section))

            # Extract interruption section
            interruption_section = match.group(2)
            sections.append(('interruption', interruption_section))

            # Continue with remainder
            remaining = remaining[match.end():]

        # Create loop wrapper functions
        def make_loop(content, is_first):
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
                # Clean path: partial block
                if is_first:
                    # First loop declares _xi
                    return f'''
                int _xi = 0;
                #pragma omp simd
                for (int _xj = 0; _xj < _nx % BLK_SZ; _xj++)
                {{
                    {content}
                }}'''
                else:
                    # Subsequent loops reuse _xi
                    return f'''
                #pragma omp simd
                for (int _xj = 0; _xj < _nx % BLK_SZ; _xj++)
                {{
                    {content}
                }}'''

        # Assemble sections: prep loops and interruptions alternate
        # prep -> interruption -> prep -> interruption -> ... -> prep
        result = []
        first_prep = True
        for section_type, content in sections:
            if section_type == 'prep':
                result.append(make_loop(content, first_prep))
                first_prep = False
            else:  # interruption
                # Interruptions go OUTSIDE loops, just append directly
                result.append(content)

        return '\n'.join(result)

    # ========== IKP Transformation Pipeline ==========

    def _render_body_preamble(self, body):
        """
        Override to hoist declarations for IKP mode.

        Calls base class for standard dereferencing, then performs smart
        variable analysis to only elevate variables that cross interruptions.
        """
        # Standard dereferencing (base class)
        body, preamble = super()._render_body_preamble(body)

        if self.ikp:
            # Split body into sections
            sections = self._ikp_split_into_sections(body)

            # Analyze which variables are used in which sections
            variables = self._ikp_analyze_variable_usage(sections)

            # Determine which variables need elevation (used in multiple sections)
            vars_to_elevate = {}
            vars_to_keep = {}

            for var_name, var_info in variables.items():
                # A variable needs elevation if it's used in:
                # 1. Multiple sections (prep OR interruption)
                # 2. OR if it's used in any interruption (interruptions are outside loops)

                sections_used = var_info['used_in']
                interruption_sections = [idx for idx in sections_used
                                        if idx < len(sections) and sections[idx][0] == 'interruption']

                # If used in an interruption, it MUST be elevated
                # (interruptions are outside loops, so variables must persist)
                if interruption_sections or len(sections_used) > 1:
                    vars_to_elevate[var_name] = var_info
                else:
                    vars_to_keep[var_name] = var_info

            # Track which arrays/scalars need transformation
            self._ikp_local_arrays = [name for name, info in vars_to_elevate.items()
                                     if info['type'] in ('array', 'scalar')]

            # Transform declarations for elevated variables
            declarations = []

            for var_name, var_info in vars_to_elevate.items():
                dtype = var_info['decl_info'][0]

                if var_info['type'] == 'const_array':
                    dims = var_info['decl_info'][1]
                    initializer = var_info['decl_info'][2]
                    decl = self._ikp_transform_const_decl(dtype, var_name, dims, initializer)
                    declarations.append(decl)

                elif var_info['type'] == 'array':
                    size = var_info['decl_info'][1]
                    decl = self._ikp_transform_local_decl(dtype, var_name, size)
                    declarations.append(decl)

                elif var_info['type'] == 'scalar':
                    # Elevate scalar to array[BLK_SZ]
                    decl = f'    {dtype} {var_name}[BLK_SZ];'
                    declarations.append(decl)

            # Remove elevated variable declarations from body
            for var_name, var_info in vars_to_elevate.items():
                # Remove array declarations
                body = re.sub(r'\s*fpdtype_t\s+' + re.escape(var_name) + r'\[[^\]]+\];', '', body)

                # Handle scalar declarations specially to preserve assignments
                if var_info['type'] == 'scalar':
                    # For scalars with initializers, preserve the assignment part
                    # Transform: fpdtype_t testScalar = expr; -> testScalar = expr;
                    scalar_decl_pattern = r'(\s*)fpdtype_t\s+' + re.escape(var_name) + r'\s*=\s*([^;]+);'
                    body = re.sub(scalar_decl_pattern, r'\1' + var_name + r' = \2;', body)
                    # Also remove simple declarations without initializers
                    body = re.sub(r'\s*fpdtype_t\s+' + re.escape(var_name) + r'\s*;', '', body)

                # Remove const array declarations
                body = re.sub(r'\s*const\s+fpdtype_t\s+' + re.escape(var_name) + r'\[[^\]]+\](?:\[[^\]]+\])*\s*=\s*\{[^}]*\}\s*;',
                             '', body, flags=re.DOTALL)

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
