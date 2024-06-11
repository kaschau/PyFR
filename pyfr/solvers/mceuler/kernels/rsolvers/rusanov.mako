<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.mixture_state'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>

    fpdtype_t q[${ns+2}];
    // Compute left thermodynamic properties
    fpdtype_t qhl[${3}];
    ${pyfr.expand('mixture_state', 'ul', 'q', 'qhl')};

    // Compute left fluxes + velocities
    fpdtype_t fl[${ndims}][${nvars}], vl[${ndims}];
    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'q[0]', 'vl')};

    // Right Mixture state
    fpdtype_t qhr[${3}];
    ${pyfr.expand('mixture_state', 'ur', 'q', 'qhr')};

    // Compute right fluxes + velocities
    fpdtype_t fr[${ndims}][${nvars}], vr[${ndims}];
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'q[0]', 'vr')};

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed / 2
    fpdtype_t a = 0.25*(qhl[2]+qhr[2]) + 0.5*fabs(nv);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
