<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, ql, qr, qhl, qhr, n, nf'>

    // Compute left fluxes
    fpdtype_t fl[${ndims}][${nvars}];
    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'ql')};

    // Compute right fluxes
    fpdtype_t fr[${ndims}][${nvars}];
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'qr')};

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'ql[{i}] + qr[{i}]', i=(1,ndims+1))};

    // Estimate the maximum wave speed / 2
    fpdtype_t a = 0.25*(qhl[2]+qhr[2]) + 0.5*fabs(nv);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
