<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## Point-implicit chemistry preconditioner application:
## y[upt] = out_scale * (Minv[upt] @ (in_scale * x[upt]))
<%pyfr:kernel name='ptchemapply' ndim='1'
              x='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              minv='in fpdtype_t[${str(nupts)}][${str(nvars*nvars)}]'
              y='out fpdtype_t[${str(nupts)}][${str(nvars)}]'>
% if in_scale:
    const fpdtype_t _isc[] = ${pyfr.carray(in_scale)};
% endif
% if out_scale:
    const fpdtype_t _osc[] = ${pyfr.carray(out_scale)};
% endif
    for (int _upt = 0; _upt < ${nupts}; _upt++)
    {
        fpdtype_t xs[${nvars}];
        for (int _c = 0; _c < ${nvars}; _c++)
% if in_scale:
            xs[_c] = _isc[_c]*x[_upt][_c];
% else:
            xs[_c] = x[_upt][_c];
% endif

        for (int _r = 0; _r < ${nvars}; _r++)
        {
            fpdtype_t _acc = 0.0;
            for (int _c = 0; _c < ${nvars}; _c++)
                _acc += minv[_upt][_r*${nvars} + _c]*xs[_c];
% if out_scale:
            y[_upt][_r] = _osc[_r]*_acc;
% else:
            y[_upt][_r] = _acc;
% endif
        }
    }
</%pyfr:kernel>
