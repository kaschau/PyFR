<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## CUDA IKP Macros
## Thread-cooperative operations for inner-kernel parallelism

<%pyfr:macro name='block-break' params='Break'>

 // Empty Block Break

</%pyfr:macro>


<%pyfr:macro name='loadv' params='dst, src, py:indices'>
<%
  # Cooperative vector load from 2D array into 1D array using arithmetic indexing
  # dst: destination array name (e.g., 'ucol')
  # src: source 2D array name (e.g., 'u')
  # indices: List of (row, col) tuples specifying which elements to load
  #          e.g., [(i, svar) for i in range(nupts)] for a column
  #          e.g., [(i, i) for i in range(n)] for diagonal

  nthreads = _kernel_generator.ikpnthrds
  indices_list = list(indices)
  nelem = len(indices_list)

  # Build index mapping arrays that we'll embed as constants
  # For each element, we need to know: what row and col to access
  rows = [row for row, col in indices_list]
  cols = [col for row, col in indices_list]
%>
  // Cooperative load: ${nelem} elements from 2D array (${nthreads} threads)
  {
    // Index maps: which (row,col) corresponds to each dst element
    const int _rows[${nelem}] = {${', '.join(map(str, rows))}};
    const int _cols[${nelem}] = {${', '.join(map(str, cols))}};

    int _tid_ = threadIdx.x % ${nthreads};
    for (int _i = _tid_; _i < ${nelem}; _i += ${nthreads})
    {
        dst[_i] = src[_rows[_i]][_cols[_i]];
    }
  }
</%pyfr:macro>


<%pyfr:macro name='gemv' params='c, b, py:A'>
<%
  import numpy as np

  # Matrix dimensions
  m, n = A.shape

  # Generate unique name for this matrix constant
  name = f'_ikp_mat_{id(A)}'

  # Flatten matrix to 1D array for embedding
  # Row-major layout: A[i][j] = flat_data[i*n + j]
  flat_data = A.flatten()

  # Get threads per element from generator configuration
  nthreads = _kernel_generator.ikpnthrds
%>
  ## GEMV: c = A @ b (thread-cooperative with ${nthreads} threads per element)
  ## Matrix ${m}x${n}, embedded as compile-time constant

  // Embed matrix data as compile-time constant
  __constant__ const fpdtype_t ${name}[${m * n}] = {
% for i, val in enumerate(flat_data):
      ${val}${',' if i < len(flat_data)-1 else ''}
% endfor
  };

  // Cooperatively load input vector into shared memory (if needed)
  // Note: Generator transforms 'b' references to shared memory access
  __syncthreads();

  // Each thread computes subset of output rows
  // Use threadIdx.x % nthreads to get position within cooperative group
  int _tid = threadIdx.x % ${nthreads};
  for (int _row = _tid; _row < ${m}; _row += ${nthreads}) {
      fpdtype_t _sum = 0.0;
      #pragma unroll
      for (int _col = 0; _col < ${n}; _col++) {
          _sum += ${name}[_row * ${n} + _col] * b[_col];
      }
      c[_row] = _sum;
  }
  __syncthreads();

</%pyfr:macro>
