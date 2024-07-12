<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.transform'/>
<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='rsolve' params='ul, ur, ql, qr, qhl, qhr, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];

    utl[0] = ul[0];
    utr[0] = ur[0];
    utl[${ndims + 1}] = ul[${ndims + 1}];
    utr[${ndims + 1}] = ur[${ndims + 1}];
    % for i in range(ns-1):
    utl[${Yix+n}] = ul[${Yix+n}];
    utr[${Yix+n}] = ur[${Yix+n}];
    % endfor

    ${pyfr.expand('transform_to', 'n', 'ul', 'utl', off=1)};
    % for i in range(ndims):
    ql[${i+1}] = utl[${i+1}]/utl[0];
    % endfor
    ${pyfr.expand('transform_to', 'n', 'ur', 'utr', off=1)};
    % for i in range(ndims):
    qr[${i+1}] = utr[${i+1}]/utr[0];
    % endfor

    ${pyfr.expand('rsolve_1d', 'utl', 'utr', 'ql', 'qr', 'qhl', 'qhr', 'ntf')};

    // Reset velocities
    % for i in range(ndims):
    ql[${i+1}] = ul[${i+1}]/ul[0];
    qr[${i+1}] = ur[${i+1}]/ur[0];
    % endfor

    nf[0] = ntf[0];
    nf[${ndims + 1}] = ntf[${ndims + 1}];
    % for i in range(ns-1):
    nf[${Yix+n}] = ntf[${Yix+n}];
    % endfor
    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', off=1)};
</%pyfr:macro>
