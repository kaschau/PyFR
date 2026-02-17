<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.multicomp.makoutil' name='mc'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-prims'/>

<% ns, vix, Eix, rhoix, pix, Tix = mc.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    // set right side primitive s
% for n in range(ns):
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