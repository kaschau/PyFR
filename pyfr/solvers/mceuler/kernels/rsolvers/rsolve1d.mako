<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\

<%pyfr:macro name='rsolve' params='ul, ur, ql, qr, qhl, qhr, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];
    fpdtype_t qtl[${nvars+1}], qtr[${nvars+1}];

    utl[0] = ul[0];
    qtl[0] = ql[0];
    utr[0] = ur[0];
    qtr[0] = qr[0];
    utl[${ndims + 1}] = ul[${ndims + 1}];
    qtl[${ndims + 1}] = ql[${ndims + 1}];
    utr[${ndims + 1}] = ur[${ndims + 1}];
    qtr[${ndims + 1}] = qr[${ndims + 1}];

    % for n in range(ns - 1):
      utl[${Yix + n}] = ul[${Yix + n}];
      qtl[${Yix + n}] = ql[${Yix + n}];
      utr[${Yix + n}] = ur[${Yix + n}];
      qtr[${Yix + n}] = qr[${Yix + n}];
    % endfor
    qtl[${nvars}] = ql[${nvars}];
    qtr[${nvars}] = qr[${nvars}];

    ${pyfr.expand('transform_to', 'n', 'ul', 'utl', off=1)};
    ${pyfr.expand('transform_to', 'n', 'ur', 'utr', off=1)};
    % for i in range(ndims):
      qtl[${i+1}] = utl[${i+1}]/utl[0];
      qtr[${i+1}] = utr[${i+1}]/utr[0];
    % endfor

    ${pyfr.expand('rsolve_1d', 'utl', 'utr', 'qtl', 'qtr', 'qhl', 'qhr', 'n', 'ntf')};

    nf[0] = ntf[0];
    nf[${ndims + 1}] = ntf[${ndims + 1}];
    % for n in range(ns - 1):
      nf[${Yix + n}] = ntf[${Yix + n}];
    % endfor

    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', off=1)};
</%pyfr:macro>
