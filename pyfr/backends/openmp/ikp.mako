<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## GEMV operation: c[m] = A[m][k] @ b[k]
## Kernel developer writes per-element, backend batches automatically
<%pyfr:macro name='gemv' params='c, A, b, py:m, py:k'>
// PYFR_IKP_MARKER: Enable inner-kernel parallelism / cache-blocking
// PYFR_IKP_PHASE_BOUNDARY: prep_end
// PYFR_IKP_OPERATION: gemv
// PYFR_IKP_SCRATCH_IN: b[${m}]
// PYFR_IKP_SCRATCH_OUT: c[${m}]
// PYFR_IKP_MATRIX: A[${m}][${k}]
// GEMV: c = A @ b (batched libxsmm call on all BLK_SZ elements)
// Cast integer addresses back to function pointers and call libxsmm
{
    // Function signature: void (*exec)(void *handle, const fpdtype_t *b, fpdtype_t *c)
    typedef void (*xsmm_func_t)(void*, const fpdtype_t*, fpdtype_t*);
    xsmm_func_t xsmm_exec = (xsmm_func_t)(size_t)${xsmm_exec_ptr};
    void *xsmm_handle = (void*)(size_t)${xsmm_blkptr};

    // Call libxsmm on the batched arrays
    // b[BLK_SZ][${k}], c[BLK_SZ][${m}]
    // libxsmm expects transposed layout: b^T[${k}][BLK_SZ], c^T[${m}][BLK_SZ]
    xsmm_exec(xsmm_handle, (const fpdtype_t*)b, (fpdtype_t*)c);
}
// PYFR_IKP_PHASE_BOUNDARY: gemm_end
</%pyfr:macro>

## GEMM operation: C[m][n] = A[m][k] @ B[k][n]
## Kernel developer writes per-element, backend batches automatically
<%pyfr:macro name='gemm' params='C, A, B, py:m, py:n, py:k'>
// PYFR_IKP_MARKER: Enable inner-kernel parallelism / cache-blocking
// GEMM: C = A @ B (per-element operation, auto-batched by backend)
% for i in range(m):
%   for j in range(n):
C[${i}][${j}] = ${''.join(f'{"+" if p > 0 else ""}A[{i}][{p}]*B[{p}][{j}]' for p in range(k))};
%   endfor
% endfor
</%pyfr:macro>
