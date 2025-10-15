<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.backends.openmp.ikputil' name='ikputil'/>

## GEMV operation: c = A @ b (batched libxsmm call)
## Kernel developer writes per-element, backend batches automatically
<%pyfr:macro name='gemv' params='c, b, py:A'>
<%
    # Call backend helper to create libxsmm kernel from the matrix data
    # A is a numpy array passed during expand()
    # Helper returns (exec_ptr, blkptr) for the JIT-compiled kernel
    # Note: Mako automatically passes context as first argument to namespace functions
    exec_ptr, blkptr = ikputil.create_gemv_kernel(A)
%>
// PYFR_IKP_MARKER: Enable inner-kernel parallelism / cache-blocking
// GEMV: c = A @ b (batched libxsmm call on all BLK_SZ elements)
{
    typedef void (*xsmm_func_t)(void*, const fpdtype_t*, fpdtype_t*);
    xsmm_func_t xsmm_exec = (xsmm_func_t)(size_t)${exec_ptr};
    void *xsmm_handle = (void*)(size_t)${blkptr};
    xsmm_exec(xsmm_handle, (const fpdtype_t*)b, (fpdtype_t*)c);
}
</%pyfr:macro>