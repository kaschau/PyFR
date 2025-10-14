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
