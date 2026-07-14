<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.chem.net-rate-jacobian'/>

## Point-implicit chemistry preconditioner build: for each solution point
## form A = I - gdt*d(src)/du from the analytic chemistry Jacobian and
## invert it in place via Gauss-Jordan with partial pivoting.  The blocks
## are nvars x nvars, stored row-major per point.
<%pyfr:kernel name='ptchemprecond' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              minv='out fpdtype_t[${str(nupts)}][${str(nvars*nvars)}]'
              gdt='scalar fpdtype_t'>
    for (int _upt = 0; _upt < ${nupts}; _upt++)
    {
        fpdtype_t uu[${nvars}];
        for (int _v = 0; _v < ${nvars}; _v++)
            uu[_v] = u[_upt][_v];

        fpdtype_t q[${nvars + 2}], qh[${4 + mcf.ns}];
        ${pyfr.expand('stateFrom-cons', 'uu', 'q', 'qh')};

        fpdtype_t jblk[${nvars*nvars}];
        ${pyfr.expand('net_rate_jacobian', 'q', 'qh', 'jblk')};

        // A = I - gdt*J; M = I
        fpdtype_t A[${nvars}][${nvars}], Mi[${nvars}][${nvars}];
        for (int _r = 0; _r < ${nvars}; _r++)
            for (int _c = 0; _c < ${nvars}; _c++)
            {
                A[_r][_c] = ((_r == _c) ? 1.0 : 0.0)
                          - gdt*jblk[_r*${nvars} + _c];
                Mi[_r][_c] = (_r == _c) ? 1.0 : 0.0;
            }

        // Gauss-Jordan with partial pivoting
        for (int _p = 0; _p < ${nvars}; _p++)
        {
            int _piv = _p;
            fpdtype_t _amax = fabs(A[_p][_p]);
            for (int _r = _p + 1; _r < ${nvars}; _r++)
                if (fabs(A[_r][_p]) > _amax)
                {
                    _amax = fabs(A[_r][_p]);
                    _piv = _r;
                }

            if (_piv != _p)
                for (int _c = 0; _c < ${nvars}; _c++)
                {
                    fpdtype_t _t = A[_p][_c];
                    A[_p][_c] = A[_piv][_c];
                    A[_piv][_c] = _t;
                    _t = Mi[_p][_c];
                    Mi[_p][_c] = Mi[_piv][_c];
                    Mi[_piv][_c] = _t;
                }

            // A singular pivot degrades to (scaled) identity behaviour
            fpdtype_t _d = 1.0/((_amax > ${fpdtype_min}) ? A[_p][_p] : 1.0);
            for (int _c = 0; _c < ${nvars}; _c++)
            {
                A[_p][_c] *= _d;
                Mi[_p][_c] *= _d;
            }

            for (int _r = 0; _r < ${nvars}; _r++)
                if (_r != _p)
                {
                    fpdtype_t _f = A[_r][_p];
                    for (int _c = 0; _c < ${nvars}; _c++)
                    {
                        A[_r][_c] -= _f*A[_p][_c];
                        Mi[_r][_c] -= _f*Mi[_p][_c];
                    }
                }
        }

        for (int _r = 0; _r < ${nvars}; _r++)
            for (int _c = 0; _c < ${nvars}; _c++)
                minv[_upt][_r*${nvars} + _c] = Mi[_r][_c];
    }
</%pyfr:kernel>
