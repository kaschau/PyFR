<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.backends.base.ikp'/>

## Metal IKP Macros
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

  indices_list = list(indices)
  nelem = len(indices_list)

  # Build index mapping arrays that we'll embed as constants
  rows = [row for row, col in indices_list]
  cols = [col for row, col in indices_list]
%>
  // Cooperative load: ${nelem} elements from 2D array (_tpitg.y threads)
  {
    // Index maps: which (row,col) corresponds to each dst element
    const int _rows[${nelem}] = {${', '.join(map(str, rows))}};
    const int _cols[${nelem}] = {${', '.join(map(str, cols))}};

    uint _tid_ = _tpitg.y;
    for (int _i = _tid_; _i < ${nelem}; _i += 8)
    {
        dst[_i] = src[_rows[_i]][_cols[_i]];
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);
</%pyfr:macro>


<%pyfr:macro name='gemv' params='c, b, py:A'>
<%
  import numpy as np

  # Matrix dimensions
  m, n = A.shape

  # Generate unique name for this matrix constant
  name = f'_ikp_mat_{id(A)}'

  # Flatten matrix to 1D array for embedding
  flat_data = A.flatten()
%>
  ## GEMV: c = A @ b (thread-cooperative with 8 threads per element)
  ## Matrix ${m}x${n}, embedded as compile-time constant

  // Embed matrix data as compile-time constant (const, not constant address space)
  const fpdtype_t ${name}[${m * n}] = {
% for i, val in enumerate(flat_data):
      ${val}${',' if i < len(flat_data)-1 else ''}
% endfor
  };

  // Each thread computes subset of output rows
  uint _tid = _tpitg.y;
  for (int _row = _tid; _row < ${m}; _row += 8) {
      fpdtype_t _sum = 0.0;
      for (int _col = 0; _col < ${n}; _col++) {
          _sum += ${name}[_row * ${n} + _col] * b[_col];
      }
      c[_row] = _sum;
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

</%pyfr:macro>


<%pyfr:macro name='square_arr' params='dst, src, py:n'>
  // Cooperatively compute dst[i] = src[i] * src[i]
  uint _tid = _tpitg.y;
  for (int _i = _tid; _i < ${n}; _i += 8) {
      dst[_i] = src[_i] * src[_i];
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

</%pyfr:macro>


<%pyfr:macro name='reduce_sum' params='_redbuf, result, src, py:n'>
<%
  import math
  ythrds = 8  # block2d[1] for Metal backend
%>
  // Cooperative reduction: result = sum(src[i] for i in range(n))
  // Works for any n (independent of block configuration)
  // Uses threadgroup memory buffer _redbuf for ${ythrds} threads per element
  uint _tid = _tpitg.y;

  // Each thread accumulates its subset
  fpdtype_t _psum = 0.0;
  for (int _i = _tid; _i < ${n}; _i += ${ythrds}) {
      _psum += src[_i];
  }
  _redbuf[_tid] = _psum;
  threadgroup_barrier(mem_flags::mem_threadgroup);

  // Tree reduction in threadgroup memory
% for step in range(int(math.log2(ythrds)), 0, -1):
<%
    stride = 2 ** (step - 1)
%>
  if (_tid < ${stride}) {
      _redbuf[_tid] += _redbuf[_tid + ${stride}];
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);
% endfor

  // All threads read the result
  result = _redbuf[0];

</%pyfr:macro>


<%pyfr:macro name='reduce_sum_masked' params='_redbuf, result, src, mask, py:n'>
<%
  import math
  ythrds = 8  # block2d[1] for Metal backend
%>
  // Cooperative masked reduction: result = sum(src[i] where mask[i])
  // Works for any n (independent of block configuration)
  // Uses threadgroup memory buffer _redbuf for ${ythrds} threads per element
  uint _tid = _tpitg.y;

  // Each thread accumulates its subset
  fpdtype_t _psum = 0.0;
  for (int _i = _tid; _i < ${n}; _i += ${ythrds}) {
      if (mask[_i]) {
          _psum += src[_i];
      }
  }
  _redbuf[_tid] = _psum;
  threadgroup_barrier(mem_flags::mem_threadgroup);

  // Tree reduction in threadgroup memory
% for step in range(int(math.log2(ythrds)), 0, -1):
<%
    stride = 2 ** (step - 1)
%>
  if (_tid < ${stride}) {
      _redbuf[_tid] += _redbuf[_tid + ${stride}];
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);
% endfor

  // All threads read the result
  result = _redbuf[0];

</%pyfr:macro>
