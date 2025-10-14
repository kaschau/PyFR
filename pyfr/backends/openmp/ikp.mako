<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## GEMV operation: c = A @ b (batched libxsmm call)
## Kernel developer writes per-element, backend batches automatically
<%pyfr:macro name='gemv' params='c, b'>
// PYFR_IKP_MARKER: Enable inner-kernel parallelism / cache-blocking
// GEMV: c = A @ b (batched libxsmm call on all BLK_SZ elements)
{
    typedef void (*xsmm_func_t)(void*, const fpdtype_t*, fpdtype_t*);
    xsmm_func_t xsmm_exec = (xsmm_func_t)(size_t)${xsmm_exec_ptr};
    void *xsmm_handle = (void*)(size_t)${xsmm_blkptr};
    xsmm_exec(xsmm_handle, (const fpdtype_t*)b, (fpdtype_t*)c);
}
</%pyfr:macro>