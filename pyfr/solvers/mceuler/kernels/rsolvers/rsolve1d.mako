<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.mixture_state'/>

<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];

    utl[0] = ul[0];
    utr[0] = ur[0];
    utl[${ndims + 1}] = ul[${ndims + 1}];
    utr[${ndims + 1}] = ur[${ndims + 1}];

    ${pyfr.expand('transform_to', 'n', 'ul', 'utl', off=1)};
    ${pyfr.expand('transform_to', 'n', 'ur', 'utr', off=1)};


    // Left Mixture state
    fpdtype_t pl, Tl, Yl[${ns}];
    fpdtype_t Rmixl, cpmixl;
    ${pyfr.expand('mixture_state', 'ul', 'pl', 'Tl', 'Yl', 'Rmixl', 'cpmixl')};

    // Right Mixture state
    fpdtype_t pr, Tr, Yr[${ns}];
    fpdtype_t Rmixr, cpmixr;
    ${pyfr.expand('mixture_state', 'ur', 'pr', 'Tr', 'Yr', 'Rmixr', 'cpmixr')};

    ${pyfr.expand('rsolve_1d', 'utl', 'utr', 'ntf')};

    nf[0] = ntf[0];
    nf[${nvars - 1}] = ntf[${nvars - 1}];
    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', off=1)};
</%pyfr:macro>
