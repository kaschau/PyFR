<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
${pyfr.eos_check('rusanov', fluid, 'mc-cpg', 'mc-tpg')}\
<%include file='pyfr.solvers.mceuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right primitive state
    ${fluid.decl('ul', 'v, p, a', suffix='l')}
    ${fluid.decl('ur', 'v, p, a', suffix='r')}

    // Compute the left and right fluxes
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    ${pyfr.expand('inviscid_flux', 'ul', 'pl', 'vl', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'pr', 'vr', 'fr')};

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${' + '.join(f'n[{i}]*(vl[{i}] + vr[{i}])'
                                for i in range(ndims))};

    // Estimate the maximum wave speed / 2
    fpdtype_t a = 0.25*(al + ar) + 0.25*fabs(nv);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
