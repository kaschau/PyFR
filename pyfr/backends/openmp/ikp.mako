<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## GEMV operation: c[m] = A[m][k] @ b[k]
## Kernel developer writes per-element, backend batches automatically
<%pyfr:macro name='gemv' params='c, A, b, py:m, py:k'>
// PYFR_IKP_MARKER: Enable inner-kernel parallelism / cache-blocking
// GEMV: c = A @ b (batched libxsmm call on all BLK_SZ elements)
{
    typedef void (*xsmm_func_t)(void*, const fpdtype_t*, fpdtype_t*);
    xsmm_func_t xsmm_exec = (xsmm_func_t)(size_t)${xsmm_exec_ptr};
    void *xsmm_handle = (void*)(size_t)${xsmm_blkptr};
    xsmm_exec(xsmm_handle, (const fpdtype_t*)b, (fpdtype_t*)c);
}
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
