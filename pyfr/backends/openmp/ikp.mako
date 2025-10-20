<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.backends.openmp.ikputil' name='ikputil'/>

<%pyfr:macro name='block-break' params='Break'>

 // Empty Block Break

</%pyfr:macro>


<%pyfr:macro name='loadv' params='dst, src, py:indices'>
<%
  # Sequential vector load from 2D array into 1D array
  # dst: destination array name (e.g., 'ucol')
  # src: source 2D array name (e.g., 'u')
  # indices: List of (row, col) tuples specifying which elements to load
  #          e.g., [(i, svar) for i in range(nupts)] for a column
  #          e.g., [(i, i) for i in range(n)] for diagonal
  indices_list = list(indices)
%>
  // Sequential load: ${len(indices_list)} elements from 2D array
% for dst_idx, (row, col) in enumerate(indices_list):
  dst[${dst_idx}] = src[${row}][${col}];
% endfor
</%pyfr:macro>


<%pyfr:macro name='gemv' params='c, b, py:A'>
<%
  # Call backend helper to create libxsmm kernel from the matrix data
  # A is a numpy array passed during expand()
  # Helper returns (exec_ptr, blkptr) for the JIT-compiled kernel
  # Note: Mako automatically passes context as first argument to namespace functions
  exec_ptr, blkptr = ikputil.create_gemv_kernel(A)
%>
  ## GEMV: c = A @ b (batched libxsmm call on all BLK_SZ elements)
  // HACK
  // NO LOOP

  typedef void (*xsmm_func_t)(void*, const fpdtype_t*, fpdtype_t*);
  xsmm_func_t xsmm_exec = (xsmm_func_t)(size_t)${exec_ptr};
  void *xsmm_handle = (void*)(size_t)${blkptr};
  xsmm_exec(xsmm_handle, (const fpdtype_t*)b, (fpdtype_t*)c);

</%pyfr:macro>