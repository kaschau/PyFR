"""
OpenMP-specific IKP utilities for Mako templates.

These functions are called during template rendering to create
libxsmm kernels for IKP interruptions.
"""

import numpy as np


def create_gemv_kernel(context, matrix):
    """
    Create a libxsmm GEMV kernel for IKP cache-blocking.

    Called from ikp.mako during template rendering to create the
    JIT-compiled libxsmm kernel.

    Parameters
    ----------
    context : Mako context
        Template rendering context (automatically passed by Mako)
    matrix_np : numpy.ndarray
        The matrix data (e.g., invvdm as numpy array)

    Returns
    -------
    exec_ptr : int
        Function pointer to libxsmm executor
    blkptr : int
        Handle to JIT-compiled libxsmm kernel
    """
    # Get backend from context
    backend = context['_backend']

    # Get matrix dimensions
    m, k = matrix.shape

    # Get block size from backend
    blk_sz = backend.csubsz

    # Create a const_matrix from the numpy array
    const_matrix = backend.const_matrix(matrix, tags={'align'})

    # Create scratch matrices for the batched GEMV operation
    # Input: b[k][blk_sz], Output: c[m][blk_sz]
    scratch_b = backend.matrix((k, blk_sz), tags={'align'})
    scratch_c = backend.matrix((m, blk_sz), tags={'align'})

    # Use the backend's xsmm kernel provider to create the mul kernel
    ikp_gemv = backend.kernel('mul', const_matrix, scratch_b, out=scratch_c)

    # Extract the function pointer and kernel handle from the kernel
    gemv_kargs = ikp_gemv.kernel.kargs
    exec_ptr = gemv_kargs.arg0  # libxsmm executor function pointer
    blkptr = gemv_kargs.arg1     # JIT-compiled kernel handle

    # Store kernel reference to prevent garbage collection
    # The kernel object must stay alive for the JIT pointers to remain valid
    backend._ikp_kernels.append(ikp_gemv)

    return exec_ptr, blkptr
