<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-prims'/>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    // set right side primitive s
% for n in range(mcf.ns):
    qr[${n}] = ql[${n}];
% endfor

% for i in range(ndims):
    qr[${i + vix}] = ql[${i + vix}];
% endfor

    // fix pressure
    qr[${pix}] = ${c['p']};

    qr[${Tix}] = ql[${Tix}];

${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>