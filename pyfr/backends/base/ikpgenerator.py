"""
Base IKP (Inner-Kernel Parallelism) generator support.

IKP enables backend-specific optimizations for heavy pointwise operations:
- OpenMP: Cache-blocking with sequential element loops + libxsmm batched calls
- GPU: Thread cooperation with shared memory + WMMA/tensor cores

This module provides common pattern matching and tracking for IKP transformations.
Backend-specific classes implement the actual transformation logic.
"""

import re


class IKPKernelGeneratorMixin:
    """
    Mixin class providing common IKP transformation utilities.

    This class handles:
    - Pattern matching for constant and local array declarations
    - Tracking which arrays need IKP transformation
    - Common regex patterns

    Backend-specific subclasses implement:
    - How to transform declarations (stack vs shared memory)
    - How to transform references (ELEM_IDX, threadIdx.x, etc.)
    - How to wrap body (sequential loop vs thread parallelism)
    """

    # Regex patterns for declaration matching (common to all backends)
    _CONST_ARRAY_PATTERN = r'\s*const\s+(fpdtype_t)\s+(\w+)(\[[^\]]+\](?:\[[^\]]+\])*)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\})\s*;'
    _LOCAL_ARRAY_PATTERN = r'\s*(fpdtype_t)\s+(\w+)\[([^\]]+)\];'
    _LOCAL_SCALAR_PATTERN = r'\s*(fpdtype_t)\s+(\w+)\s*(?:=\s*[^;]+)?;'

    def _ikp_split_into_sections(self, body):
        """
        Split IKP body into alternating prep/interruption sections.

        Returns:
            list of (section_type, content) tuples
            section_type is 'prep' or 'interruption'
        """
        sections = []
        remaining = body

        while True:
            # Find next interruption
            match = re.search(
                r'(.*?)// PYFR_IKP_INTERRUPTION_START\n(.*?)// PYFR_IKP_INTERRUPTION_END',
                remaining, flags=re.DOTALL
            )

            if not match:
                # No more interruptions, rest is final prep/proc section
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

        return sections

    def _ikp_analyze_variable_usage(self, sections):
        """
        Analyze which variables are used in which sections.

        Returns:
            dict: {var_name: {'type': 'array'/'scalar'/'const_array',
                              'declared_in': section_index,
                              'used_in': set of section_indices,
                              'decl_info': (dtype, size/dims, initializer)}}
        """
        variables = {}

        for idx, (section_type, content) in enumerate(sections):
            # Skip interruptions - they don't have local declarations
            if section_type == 'interruption':
                continue

            # Find constant array declarations
            for match in re.finditer(self._CONST_ARRAY_PATTERN, content, re.DOTALL):
                dtype, name, dims, initializer = match.groups()
                if name not in variables:
                    variables[name] = {
                        'type': 'const_array',
                        'declared_in': idx,
                        'used_in': set(),
                        'decl_info': (dtype, dims, initializer)
                    }

            # Find local array declarations
            for match in re.finditer(self._LOCAL_ARRAY_PATTERN, content):
                dtype, name, size = match.groups()
                if name not in variables:
                    variables[name] = {
                        'type': 'array',
                        'declared_in': idx,
                        'used_in': set(),
                        'decl_info': (dtype, size, None)
                    }

            # Find local scalar declarations (but filter out arrays)
            # Remove array declarations first to avoid false matches
            content_no_arrays = re.sub(self._LOCAL_ARRAY_PATTERN, '', content)
            content_no_arrays = re.sub(self._CONST_ARRAY_PATTERN, '', content_no_arrays, flags=re.DOTALL)

            for match in re.finditer(self._LOCAL_SCALAR_PATTERN, content_no_arrays):
                dtype, name = match.groups()
                # Skip common C keywords and types
                if name in {'void', 'int', 'char', 'float', 'double', 'if', 'for', 'while', 'return'}:
                    continue
                if name not in variables:
                    variables[name] = {
                        'type': 'scalar',
                        'declared_in': idx,
                        'used_in': set(),
                        'decl_info': (dtype, None, None)
                    }

        # Now find usages of each variable across all sections
        for var_name in variables.keys():
            for idx, (section_type, content) in enumerate(sections):
                # Check if variable is used in this section (look for var_name as a word)
                if re.search(r'\b' + re.escape(var_name) + r'\b', content):
                    variables[var_name]['used_in'].add(idx)

        return variables

    def _ikp_find_declarations(self, body):
        """
        Find constant and local array declarations in kernel body.

        Returns:
            tuple: (const_arrays, local_arrays)
                const_arrays: list of (dtype, name, dims, initializer)
                local_arrays: list of (dtype, name, size)
        """
        const_arrays = []
        local_arrays = []

        # Find constant array declarations
        for match in re.finditer(self._CONST_ARRAY_PATTERN, body, re.DOTALL):
            dtype, name, dims, initializer = match.groups()
            const_arrays.append((dtype, name, dims, initializer))

        # Find local array declarations
        for match in re.finditer(self._LOCAL_ARRAY_PATTERN, body):
            dtype, name, size = match.groups()
            local_arrays.append((dtype, name, size))

        return const_arrays, local_arrays

    def _ikp_remove_declarations(self, body):
        """
        Remove constant and local array declarations from body.

        These will be hoisted to preamble with backend-specific transformations.
        """
        # Remove constant arrays
        body = re.sub(self._CONST_ARRAY_PATTERN, '', body, flags=re.DOTALL)

        # Remove local arrays
        body = re.sub(self._LOCAL_ARRAY_PATTERN, '', body)

        return body

    # Backend-specific methods (must be implemented by subclasses)

    def _ikp_transform_const_decl(self, dtype, name, dims, initializer):
        """
        Transform constant array declaration for IKP.

        OpenMP: Keep as-is (shared across elements)
        GPU: Might need __constant__ qualifier
        """
        raise NotImplementedError("Backend must implement _ikp_transform_const_decl")

    def _ikp_transform_local_decl(self, dtype, name, size):
        """
        Transform local array declaration for IKP.

        OpenMP: arr[N] -> arr[BLK_SZ][N] (stack allocation)
        GPU: arr[N] -> __shared__ arr[BLOCK_DIM][N] (shared memory)
        """
        raise NotImplementedError("Backend must implement _ikp_transform_local_decl")

    def _ikp_transform_array_ref(self, arr_name, body):
        """
        Transform array references in IKP sections.

        OpenMP: arr[i] -> arr[ELEM_IDX][i]
        GPU: arr[i] -> arr[threadIdx.x][i]
        """
        raise NotImplementedError("Backend must implement _ikp_transform_array_ref")

    def _ikp_transform_kernel_args(self, body):
        """
        Transform kernel argument references in IKP sections.

        OpenMP: X_IDX_AOSOA(...) -> ELEM_IDX
        GPU: Similar but might use threadIdx.x directly
        """
        raise NotImplementedError("Backend must implement _ikp_transform_kernel_args")

    def _ikp_elem_idx_macro(self):
        """
        Return the ELEM_IDX macro definition for this backend.

        OpenMP: #define ELEM_IDX _elem
        GPU: #define ELEM_IDX threadIdx.x
        """
        raise NotImplementedError("Backend must implement _ikp_elem_idx_macro")

    def _ikp_wrap_body(self, body, nelem_expr):
        """
        Wrap transformed body in backend-specific parallelization.

        OpenMP: for (_elem = 0; _elem < nelem_expr; _elem++) { body }
        GPU: Just body (each thread is an element)
        """
        raise NotImplementedError("Backend must implement _ikp_wrap_body")
