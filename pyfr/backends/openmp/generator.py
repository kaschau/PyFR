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

    """

    # ========== IKP Mixin Implementation (OpenMP-Specific) ==========

    def _ikp_transform_local_decl(self, localvar):
        """
        Transform local array declaration for OpenMP IKP.

        Strategy: Only flatten the first dimension (element parallelism)
        - 1D: arr[N] -> arr[BLK_SZ*N]
        - 2D: arr[M][N] -> arr[BLK_SZ*M][N]
        - 3D: arr[M][N][P] -> arr[BLK_SZ*M][N][P]
        - etc.
        """
        if localvar.ncdim == 0:
            raise ValueError(f'Scalars should not call _ikp_transform_local_decl')

        # Flatten first dimension, keep rest intact
        first_dim = localvar.cdims[0]
        rest_dims = ''.join(f'[{d}]' for d in localvar.cdims[1:])

        return f'    {localvar.dtype} {localvar.name}[BLK_SZ*{first_dim}]{rest_dims};'

    def _ikp_transform_array_ref(self, localvar, body):
        """
        Transform local variable references for OpenMP IKP.

        Strategy: Only transform the first index (element parallelism)
        - Scalars: x -> x[X_IDX]
        - 1D: arr[i] -> arr[i*BLK_SZ + X_IDX]
        - 2D: arr[i][j] -> arr[i*BLK_SZ + X_IDX][j]
        - 3D: arr[i][j][k] -> arr[i*BLK_SZ + X_IDX][j][k]
        - etc.
        """
        var_name = localvar.name

        if localvar.isscalar:
            # Transform scalar -> scalar[X_IDX]
            scalar_pattern = rf'\b{var_name}\b(?!\s*\[)'
            return re.sub(scalar_pattern, rf'{var_name}[X_IDX]', body)
        else:
            # For arrays: transform first index only
            # Build pattern that matches all dimensions: arr[i][j][k]...
            # Capture groups: (i), (j), (k), ...
            bracket_pattern = r'\[([^\]]+)\]'
            pattern = rf'\b{var_name}' + bracket_pattern * localvar.ncdim

            # Build replacement: arr[(i)*BLK_SZ + X_IDX][j][k]...
            # First index gets transformed, rest stay the same
            replacement = rf'{var_name}[(\1)*BLK_SZ + X_IDX]'
            for i in range(2, localvar.ncdim + 1):
                replacement += rf'[\{i}]'

            return re.sub(pattern, replacement, body)

    def _ikp_wrap_body(self, body, nelem_expr):
        """
        Wrap IKP body in standard _xi/_xj loop structure for parallel execution.

        Splits body into alternating prep/interruption sections:
        - Prep sections (per-element code) - INSIDE element loops
        - Interruption sections (batched ops) - OUTSIDE loops
        Handles N interruptions automatically.
        """
        # Split body on interruption markers using common method
        sections = self._ikp_split_into_sections(body)

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

            for var_name, var_info in variables.items():
                sections_used = var_info['used_in']
                interruption_sections = [idx for idx in sections_used
                                        if idx < len(sections) and sections[idx][0] == 'interruption']

                # If used in an interruption, it MUST be elevated
                # (interruptions are outside loops, so variables must persist)
                if interruption_sections or len(sections_used) > 1:
                    vars_to_elevate[var_name] = var_info

            # Track which arrays/scalars need transformation (store IKPLocalVar objects)
            self._ikp_local_vars = {name: info['localvar']
                                    for name, info in vars_to_elevate.items()}

            # Transform declarations for elevated variables
            declarations = []

            for var_name, var_info in vars_to_elevate.items():
                localvar = var_info['localvar']

                if localvar.isarray:
                    # Transform array: arr[N] -> arr[BLK_SZ*N] or arr[M][N] -> arr[BLK_SZ*M][N]
                    decl = self._ikp_transform_local_decl(localvar)
                    declarations.append(decl)
                else:
                    # Elevate scalar to array[BLK_SZ]
                    decl = f'    {localvar.dtype} {var_name}[BLK_SZ];'
                    declarations.append(decl)

            # Remove elevated variable declarations from body
            for var_name, var_info in vars_to_elevate.items():
                localvar = var_info['localvar']
                dtype_escaped = re.escape(localvar.dtype)
                name_escaped = re.escape(var_name)

                if localvar.isarray:
                    # Remove array declarations (any type, any dimensions, with optional initializer)
                    # Matches: int arr[5]; fpdtype_t mat[3][4][2] = {...}; etc.
                    pattern = rf'\s*{dtype_escaped}\s+{name_escaped}(?:\[\d+\])+\s*(?:=\s*\{{[^}}]*\}})?\s*;'
                    body = re.sub(pattern, '', body)
                else:
                    # Handle scalar declarations specially to preserve assignments
                    # Transform: int x = expr; -> x = expr;
                    scalar_decl_pattern = rf'(\s*){dtype_escaped}\s+{name_escaped}\s*=\s*([^;]+);'
                    body = re.sub(scalar_decl_pattern, rf'\1{var_name} = \2;', body)
                    # Also remove simple declarations without initializers
                    body = re.sub(rf'\s*{dtype_escaped}\s+{name_escaped}\s*;', '', body)

            # Transform array references for elevated variables
            for localvar in self._ikp_local_vars.values():
                body = self._ikp_transform_array_ref(localvar, body)

            # Add transformed declarations to preamble
            if declarations:
                preamble = '\n'.join(declarations) + '\n' + preamble

        return body, preamble

    # ========== Standard OpenMP Generator Methods ==========

    def render(self):
        kargdefn, kargassn = self._render_args('args')

        if self.ikp:
            # IKP mode: wrap body in loops (refs already transformed in _render_body_preamble)
            core = self._ikp_wrap_body(self.body, 'BLK_SZ')
            clean = self._ikp_wrap_body(self.body, '(_nx % BLK_SZ)')
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
