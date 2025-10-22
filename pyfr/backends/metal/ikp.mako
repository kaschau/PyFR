<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## Metal IKP Macros
## Thread-cooperative operations for inner-kernel parallelism

<%pyfr:macro name='block-break' params='Break'>

 // Empty Block Break

</%pyfr:macro>


<%pyfr:macro name='loadv' params='dst, src, py:indices'>
<%
  indices_list = list(indices)
  nelem = len(indices_list)
  rows = [row for row, col in indices_list]
  cols = [col for row, col in indices_list]

  # Detect if rows are sequential starting from 0 and cols are all the same
  rows_sequential = (rows == list(range(nelem)))
  cols_constant = (len(set(cols)) == 1)

  nthreads_y = _kernel_generator.block2d[1]
%>
% if rows_sequential and cols_constant:
<%
  const_col = cols[0]
%>
  // Cooperative load: ${nelem} elements from 2D array (${nthreads_y} threads, shared memory)
  // Optimized: sequential rows, constant column ${const_col}
  {
    uint _tid_ = _tpitg.y;
    for (int _i = _tid_; _i < ${nelem}; _i += ${nthreads_y})
    {
        dst[_i] = src[_i][${const_col}];
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);
% else:
  // Cooperative load: ${nelem} elements from 2D array (${nthreads_y} threads, shared memory)
  {
    const int _rows[${nelem}] = {${', '.join(map(str, rows))}};
    const int _cols[${nelem}] = {${', '.join(map(str, cols))}};

    uint _tid_ = _tpitg.y;
    for (int _i = _tid_; _i < ${nelem}; _i += ${nthreads_y})
    {
        dst[_i] = src[_rows[_i]][_cols[_i]];
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);
% endif
</%pyfr:macro>


<%pyfr:macro name='gemv' params='c, b, A, py:m, py:n'>
<%
  nthreads_y = _kernel_generator.block2d[1]
%>
  // GEMV: c = A @ b (thread-cooperative with ${nthreads_y} threads per element, shared memory)

  uint _tid = _tpitg.y;
  for (int _row = _tid; _row < ${m}; _row += ${nthreads_y}) {
      fpdtype_t _sum = 0.0;
      for (int _col = 0; _col < ${n}; _col++) {
          _sum += A[_row][_col] * b[_col];
      }
      c[_row] = _sum;
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

</%pyfr:macro>


<%pyfr:macro name='gemv_simdgroup' params='c, b, A, py:m, py:n'>
#include <metal_simdgroup_matrix>

<%
  import math
  ntiles_row = int(math.ceil(m / 8.0))
  ntiles_col = int(math.ceil(n / 8.0))
  n_padded = ntiles_col * 8
  m_padded = ntiles_row * 8
%>

  threadgroup fpdtype_t _b_matrix[${n_padded * 8}];
  {
    uint _tid = _tpitg.y;
    for (int _i = _tid; _i < ${n_padded * 8}; _i += 32) {
      int _row = _i / 8;
      int _col = _i % 8;
      _b_matrix[_i] = (_col == 0 && _row < ${n}) ? b[_row] : 0.0;
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

  threadgroup fpdtype_t _c_matrix[${m_padded * 8}];
  {
    for (int _tile_row = 0; _tile_row < ${ntiles_row}; _tile_row++) {
      simdgroup_float8x8 _acc = make_filled_simdgroup_matrix<float, 8, 8>(0.0);

      for (int _tile_k = 0; _tile_k < ${ntiles_col}; _tile_k++) {
        simdgroup_float8x8 _A_tile, _B_tile;

        simdgroup_load(_A_tile,
                      &A[_tile_row*8][_tile_k*8],
                      ${n},
                      ulong2(0, 0),
                      false);

        simdgroup_load(_B_tile,
                      &_b_matrix[_tile_k*${8*8}],
                      8,
                      ulong2(0, 0),
                      false);

        simdgroup_multiply_accumulate(_acc, _A_tile, _B_tile, _acc);
      }

      simdgroup_store(_acc,
                     &_c_matrix[_tile_row*${8*8}],
                     8,
                     ulong2(0, 0),
                     false);
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

  {
    uint _tid = _tpitg.y;
    for (int _i = _tid; _i < ${m}; _i += 32) {
      c[_i] = _c_matrix[_i * 8];
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

</%pyfr:macro>


<%pyfr:macro name='square_arr' params='dst, src, py:n'>
<%
  nthreads_y = _kernel_generator.block2d[1]
%>
  // Cooperatively compute dst[i] = src[i] * src[i] (shared memory)
  uint _tid = _tpitg.y;
  for (int _i = _tid; _i < ${n}; _i += ${nthreads_y}) {
      dst[_i] = src[_i] * src[_i];
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

</%pyfr:macro>


<%pyfr:macro name='reduce_sum' params='_redbuf, result, src, py:n'>
<%
  import math
  ythrds = _kernel_generator.block2d[1]
%>
  // Cooperative reduction using threadgroup memory
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
  ythrds = _kernel_generator.block2d[1]
%>
  // Cooperative masked reduction using threadgroup memory
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
