<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.mixture_state'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>

    // Left Mixture state
    fpdtype_t pl, Tl, Yl[${ns}];
    fpdtype_t Rmixl, cpmixl;
    ${pyfr.expand('mixture_state', 'ul', 'pl', 'Tl', 'Yl', 'Rmixl', 'cpmixl')};

    // Compute left fluxes + velocities
    fpdtype_t fl[${ndims}][${nvars}], vl[${ndims}];
    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};

    // Right Mixture state
    fpdtype_t pr, Tr, Yr[${ns}];
    fpdtype_t Rmixr, cpmixr;
    ${pyfr.expand('mixture_state', 'ur', 'pr', 'Tr', 'Yr', 'Rmixr', 'cpmixr')};

    // Compute right fluxes + velocities
    fpdtype_t fr[${ndims}][${nvars}], vr[${ndims}];
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed / 2
    fpdtype_t gammal = cpmixl/(cpmixl-Rmixl);
    fpdtype_t gammar = cpmixr/(cpmixr-Rmixr);
    fpdtype_t a = 0.25*(sqrt(gammal*Rmixl*Tl)+sqrt(gammar*Rmixr*Tr))+ 0.25*fabs(nv);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
