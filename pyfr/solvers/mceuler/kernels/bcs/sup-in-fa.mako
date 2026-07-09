<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t Tb = ${c['T']};
    fpdtype_t pb = ${c['p']};
    fpdtype_t Yb[${ns}];
% for n, spn in enumerate(sp_names):
    Yb[${n}] = ${c[spn]};
% endfor

    fpdtype_t Rb = ${' + '.join(f'Yb[{n}]*{r}'
                                for n, r in enumerate(sp_R))};
    fpdtype_t hb = ${' + '.join(f"Yb[{n}]*({fluid.h_of_T(n, 'Tb')})"
                                for n in range(ns))};
    fpdtype_t rhob = pb/(Rb*Tb);

    fpdtype_t vb[${ndims}];
% for i, v in enumerate('uvw'[:ndims]):
    vb[${i}] = ${c[v]};
% endfor

% for n in range(ns):
    ur[${n}] = rhob*Yb[${n}];
% endfor
% for i in range(ndims):
    ur[${ns + i}] = rhob*vb[${i}];
% endfor
    ur[${nvars - 1}] = rhob*hb - pb
        + 0.5*rhob*(${' + '.join(f'vb[{i}]*vb[{i}]'
                                 for i in range(ndims))});
</%pyfr:macro>
