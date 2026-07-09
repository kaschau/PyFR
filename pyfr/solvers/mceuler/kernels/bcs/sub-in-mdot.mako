<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ${fluid.decl('ul', 'rho, invrho, p', suffix='wm')}

    fpdtype_t Tb = ${c['T']};
    fpdtype_t Yb[${ns}];
% for n, spn in enumerate(sp_names):
    Yb[${n}] = ${c[spn]};
% endfor

    fpdtype_t Rb = ${' + '.join(f'Yb[{n}]*{r}'
                                for n, r in enumerate(sp_R))};
    fpdtype_t rhob = pwm/(Rb*Tb);
    fpdtype_t hb = ${' + '.join(f"Yb[{n}]*({fluid.h_of_T(n, 'Tb')})"
                                for n in range(ns))};

% for n in range(ns):
    ur[${n}] = rhob*Yb[${n}];
% endfor

    // Set the momentum to achieve the requested mass flow
% for i in range(ndims):
    ur[${ns + i}] = -2*nl[${i}]*(${c['mdot-per-area']}) - ul[${ns + i}];
% endfor

    ur[${nvars - 1}] = rhob*hb - pwm
        + (0.5/rhob)*${pyfr.dot('ur[{i}]', i=(ns, ns + ndims))};
</%pyfr:macro>
