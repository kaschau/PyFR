<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, p, v, f'>
    fpdtype_t rhoE = s[${nvars - 1}];

    // Species fluxes
% for i, n in pyfr.ndrange(ndims, ns):
    f[${i}][${n}] = v[${i}]*s[${n}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${ns + j}] = s[${ns + i}]*v[${j}]${' + p' if i == j else ''};
% endfor

    // Energy flux
% for i in range(ndims):
    f[${i}][${nvars - 1}] = (rhoE + p)*v[${i}];
% endfor
</%pyfr:macro>
