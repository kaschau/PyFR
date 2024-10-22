<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\

<%pyfr:macro name='rsolve' params='ul, ur, ql, qr, qhl, qhr, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];
    fpdtype_t qtransl[${nvars + 2}], qtransr[${nvars + 2}];

    % for n in range(ns):
      utl[${n}] = ul[${n}];
      qtransl[${n}] = ql[${n}];
      utr[${n}] = ur[${n}];
      qtransr[${n}] = qr[${n}];
    % endfor

    utl[${Eix}] = ul[${Eix}];
    utr[${Eix}] = ur[${Eix}];

    qtransl[${rhoix}] = ql[${rhoix}];
    qtransr[${rhoix}] = qr[${rhoix}];
    qtransl[${pix}] = ql[${pix}];
    qtransr[${pix}] = qr[${pix}];
    qtransl[${Tix}] = ql[${Tix}];
    qtransr[${Tix}] = qr[${Tix}];


    ${pyfr.expand('transform_to', 'n', 'ul', 'utl', off=vix)};
    ${pyfr.expand('transform_to', 'n', 'ur', 'utr', off=vix)};

    % for i in range(ndims):
      qtransl[${i+vix}] = utl[${i+vix}]/qtransl[${rhoix}];
      qtransr[${i+vix}] = utr[${i+vix}]/qtransr[${rhoix}];
    % endfor

    ${pyfr.expand('rsolve_1d', 'utl', 'utr', 'qtransl', 'qtransr', 'qhl', 'qhr', 'n', 'ntf')};

    % for n in range(ns):
      nf[${n}] = ntf[${n}];
    % endfor
    nf[${Eix}] = ntf[${Eix}];

    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', off=vix)};
</%pyfr:macro>
