<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.multicomp.makoutil' name='mc'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-prims'/>

<% ns, vix, Eix, rhoix, pix, Tix = mc.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    // set right side primitive s
%  for n,spn in enumerate(c['names']):
    qr[${n}] = ${c[spn]};
%  endfor

% for i in range(ndims):
    qr[${i + vix}] = 0.0;
% endfor

    qr[${pix}] = ql[${pix}];
    qr[${Tix}] = ${c['T']};

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

    // We now have a valid density value, set momentums to reach input mdot
    fpdtype_t invrho = 1.0/qr[${rhoix}];
% for i in range(ndims):
    ur[${i + vix}] = -2*nl[${i}]*${c['mdot-per-area']} - ul[${i + vix}];
    qr[${i + vix}] = ur[${i + vix}]*invrho;
% endfor

   // We have added ke, so we need to update total energy
    ur[${Eix}] += 0.5*invrho*${pyfr.dot('ur[{i}]', i=(vix,vix + ndims))};

</%pyfr:macro>