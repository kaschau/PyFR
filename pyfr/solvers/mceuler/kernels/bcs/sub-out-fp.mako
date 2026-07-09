<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ${fluid.decl('ul', 'rho, invrho, v, Y, T, R', suffix='wo')}

    fpdtype_t pb = ${c['p']};
    fpdtype_t rhob = pb/(Rwo*Two);
    fpdtype_t hb = ${' + '.join(f"Ywo[{n}]*({fluid.h_of_T(n, 'Two')})"
                                for n in range(ns))};

% for n in range(ns):
    ur[${n}] = rhob*Ywo[${n}];
% endfor
% for i in range(ndims):
    ur[${ns + i}] = rhob*vwo[${i}];
% endfor
    ur[${nvars - 1}] = rhob*hb - pb
        + 0.5*rhob*(${' + '.join(f'vwo[{i}]*vwo[{i}]'
                                 for i in range(ndims))});
</%pyfr:macro>
