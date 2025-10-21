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

  indices_list = list(indices)
  nelem = len(indices_list)

  # Build index mapping arrays that we'll embed as constants
  # For each element, we need to know: what row and col to access
  rows = [row for row, col in indices_list]
  cols = [col for row, col in indices_list]
%>
  // Cooperative load: ${nelem} elements from 2D array (blockDim.y threads)
  {
    // Index maps: which (row,col) corresponds to each dst element
    const int _rows[${nelem}] = {${', '.join(map(str, rows))}};
    const int _cols[${nelem}] = {${', '.join(map(str, cols))}};

    int _tid_ = threadIdx.y;
    for (int _i = _tid_; _i < ${nelem}; _i += blockDim.y)
    {
        dst[_i] = src[_rows[_i]][_cols[_i]];
    }
  }
  __syncthreads();
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
%>
  ## GEMV: c = A @ b (thread-cooperative with blockDim.y threads per element)
  ## Matrix ${m}x${n}, embedded as compile-time constant

  // Embed matrix data as compile-time constant
  __constant__ const fpdtype_t ${name}[${m * n}] = {
% for i, val in enumerate(flat_data):
      ${val}${',' if i < len(flat_data)-1 else ''}
% endfor
  };

  // Each thread computes subset of output rows
  int _tid = threadIdx.y;
  for (int _row = _tid; _row < ${m}; _row += blockDim.y) {
      fpdtype_t _sum = 0.0;
      #pragma unroll
      for (int _col = 0; _col < ${n}; _col++) {
          _sum += ${name}[_row * ${n} + _col] * b[_col];
      }
      c[_row] = _sum;
  }
  __syncthreads();

</%pyfr:macro>


<%pyfr:macro name='square_arr' params='dst, src, py:n'>
  // Cooperatively compute dst[i] = src[i] * src[i]
  int _tid = threadIdx.y;
  for (int _i = _tid; _i < ${n}; _i += blockDim.y) {
      dst[_i] = src[_i] * src[_i];
  }
  __syncthreads();

</%pyfr:macro>


<%pyfr:macro name='reduce_sum' params='_redbuf, result, src, py:n'>
<%
  import math
  ythrds = 8  # block2d[1] for CUDA backend
%>
  // Cooperative reduction: result = sum(src[i] for i in range(n))
  // Works for any n (independent of block2d configuration)
  // Uses shared memory buffer _redbuf for blockDim.y=${ythrds} threads per element
  int _tid = threadIdx.y;

  // Each thread accumulates its subset
  fpdtype_t _psum = 0.0;
  for (int _i = _tid; _i < ${n}; _i += ${ythrds}) {
      _psum += src[_i];
  }
  _redbuf[_tid] = _psum;
  __syncthreads();

  // Tree reduction in shared memory
% for step in range(int(math.log2(ythrds)), 0, -1):
<%
    stride = 2 ** (step - 1)
%>
  if (_tid < ${stride}) {
      _redbuf[_tid] += _redbuf[_tid + ${stride}];
  }
  __syncthreads();
% endfor

  // All threads read the result
  result = _redbuf[0];

</%pyfr:macro>


<%pyfr:macro name='reduce_sum_masked' params='_redbuf, result, src, mask, py:n'>
<%
  import math
  ythrds = 8  # block2d[1] for CUDA backend
%>
  // Cooperative masked reduction: result = sum(src[i] where mask[i])
  // Works for any n (independent of block2d configuration)
  // Uses shared memory buffer _redbuf for blockDim.y=${ythrds} threads per element
  int _tid = threadIdx.y;

  // Each thread accumulates its subset
  fpdtype_t _psum = 0.0;
  for (int _i = _tid; _i < ${n}; _i += ${ythrds}) {
      if (mask[_i]) {
          _psum += src[_i];
      }
  }
  _redbuf[_tid] = _psum;
  __syncthreads();

  // Tree reduction in shared memory
% for step in range(int(math.log2(ythrds)), 0, -1):
<%
    stride = 2 ** (step - 1)
%>
  if (_tid < ${stride}) {
      _redbuf[_tid] += _redbuf[_tid + ${stride}];
  }
  __syncthreads();
% endfor

  // All threads read the result
  result = _redbuf[0];

</%pyfr:macro>
