<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-prims'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    fpdtype_t qr[${nvars + 2}];

%  for n,spn in enumerate(c['names']):
    qr[${n}] = ${c[spn]};
%  endfor

% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + vix}] = ${c[v]};
% endfor

    qr[${pix}] = ${c['p']};
    qr[${Tix}] = ${c['T']};

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>
