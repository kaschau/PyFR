from math import prod
import re

from pyfr.backends.base.generator import BaseKernelGenerator
from pyfr.backends.base.ikpgenerator import IKPKernelGeneratorMixin


class OpenMPIKPKernelGeneratorMixin(IKPKernelGeneratorMixin):
    """
    OpenMP-specific IKP transformation logic.

    Handles variable elevation, declaration transformation, and reference transformation
    for OpenMP's cache-blocking strategy.
    """

    def _decl_hoists(self, hvars):
        """
        Generate C declarations for all hoisted variables.

        Strategy: Flatten first dimension for block-level parallelism
        - Scalar: x -> x[BLK_SZ]
        - 1D array: arr[N] -> arr[BLK_SZ*N]
        - 2D array: arr[M][N] -> arr[BLK_SZ*M][N]
        - ND array: arr[...] -> arr[BLK_SZ*first_dim][remaining dims]

        Returns:
            list: C declaration strings
        """
        decls = []
        for hvar in hvars:
            if hvar.isscalar:
                decls.append(f'    {hvar.dtype} {hvar.name}[BLK_SZ];')
            else:
                ldim = hvar.cdims[0]
                tdims = ''.join(f'[{d}]' for d in hvar.cdims[1:])
                decls.append(f'    {hvar.dtype} {hvar.name}[BLK_SZ*{ldim}]{tdims};')

        return decls

    def _deref_hoists(self, hvars, body):
        """
        Transform references to all hoisted variables in body.

        Strategy: Only transform the first index (element parallelism)
        - Scalars: x -> x[X_IDX]
        - 1D: arr[i] -> arr[i*BLK_SZ + X_IDX]
        - 2D: arr[i][j] -> arr[i*BLK_SZ + X_IDX][j]
        - 3D: arr[i][j][k] -> arr[i*BLK_SZ + X_IDX][j][k]

        Returns:
            str: Body with transformed references
        """
        for hvar in hvars:
            if hvar.isscalar:
                pattern = rf'\b{hvar.name}\b(?!\s*\[)'
                body = re.sub(pattern, rf'{hvar.name}[X_IDX]', body)
            else:
                # Build pattern that matches all dimensions: arr[i][j][k]...
                pattern = rf'\b{hvar.name}' + r'\[([^\]]+)\]' * hvar.ncdim

                # Build replacement: arr[(i)*BLK_SZ + X_IDX][j][k]...
                replacement = rf'{hvar.name}[(\1)*BLK_SZ + X_IDX]'
                for i in range(2, hvar.ncdim + 1):
                    replacement += rf'[\{i}]'

                body = re.sub(pattern, replacement, body)

        return body

    def _ikp_wrap_body(self, body, nelem_expr):
        """
        Wrap IKP body in _xi/_xj loop structure for parallel execution.

        Sections marked 'looped' get per-element loops.
        Sections marked 'noloop' are placed outside loops (batched operations).
        """
        sections = self._split_ikp_sections(body)

        # Create loop wrapper function
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

        # Assemble sections: wrap looped sections, append noloop sections directly
        result = []
        first_loop = True
        for stype, content in sections:
            if stype == 'looped':
                result.append(make_loop(content, first_loop))
                first_loop = False
            else:  # noloop
                result.append(content)

        return '\n'.join(result)

    def _remove_hoisted_decls(self, body, hvars):
        """
        Remove hoisted variable declarations from body.

        - Arrays: Remove entire declaration
        - Scalars: Convert "int x = expr;" to "x = expr;" or remove if no initializer

        Returns:
            str: Body with declarations removed
        """
        for hvar in hvars:
            if hvar.isarray:
                # Remove array declarations (any dimensions, optional initializer)
                pattern = rf'\s*{hvar.dtype}\s+{hvar.name}(?:\[\d+\])+\s*(?:=\s*\{{[^}}]*\}})?\s*;'
                body = re.sub(pattern, '', body)
            else:
                # Preserve initializer: int x = expr; -> x = expr;
                pattern = rf'(\s*){hvar.dtype}\s+{hvar.name}\s*=\s*([^;]+);'
                body = re.sub(pattern, rf'\1{hvar.name} = \2;', body)
                # Remove declarations without initializers
                body = re.sub(rf'\s*{hvar.dtype}\s+{hvar.name}\s*;', '', body)

        return body

    def _ikp_render_body_preamble(self, body, preamble):
        """
        Main IKP transformation pipeline for OpenMP backend.

        Orchestrates the complete transformation process:
        1. Analyze variable usage across sections
        2. Determine which variables need elevation
        3. Generate elevated declarations
        4. Remove local declarations from body
        5. Transform variable references

        Returns:
            tuple: (transformed_body, transformed_preamble)
        """
        # Split body into sections and find variables needing hoisting
        secs = self._split_ikp_sections(body)
        hvars = self._find_hvars(secs)

        # Generate preamble level declarations
        predecl = self._decl_hoists(hvars)

        # Remove original declarations from body
        body = self._remove_hoisted_decls(body, hvars)

        # Transform variable references
        body = self._deref_hoists(hvars, body)

        # Add declarations to preamble
        if predecl:
            preamble = '\n'.join(predecl) + '\n' + preamble

        return body, preamble


class OpenMPKernelGenerator(OpenMPIKPKernelGeneratorMixin, BaseKernelGenerator):
    def _render_body_preamble(self, body):
        """
        Transform kernel body and generate variable declarations.

        Calls base class for standard dereferencing, then performs IKP
        transformations if enabled.
        """
        # Standard dereferencing (base class)
        body, preamble = super()._render_body_preamble(body)

        # Apply IKP transformations if enabled
        if self.ikp:
            body, preamble = self._ikp_render_body_preamble(body, preamble)

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

        return f'''
            struct {self.name}_kargs {{ {kargdefn}; }};
            void {self.name}(ixdtype_t _ib,
                             const struct {self.name}_kargs *args,
                             int _disp_mask)
            {{
                {kargassn};
                #define X_IDX (_xi + _xj)
                #define X_IDX_AOSOA(v, nv)\
                    ((_xi/SOA_SZ*(nv) + (v))*SOA_SZ + _xj)
                #define BCAST_BLK(r, c, ld) ((c) % (ld) + ((c) / (ld))*(ld)*r)
                {self.preamble}
                if (_nx - _ib*BLK_SZ >= BLK_SZ)
                {{
                    {core}
                }}
                else
                {{
                    {clean}
                }}
                #undef X_IDX
                #undef X_IDX_AOSOA
                #undef BCAST_BLK
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
