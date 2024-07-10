<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
<% ns = c['ns'] %>
<% Yix = ndims+2 %>

    // set right side primatives
    qr[0] = ql[0];
% for i in range(ndims):
    qr[${i + 1}] = 0.0;
% endfor
    qr[${ndims+1}] = ${c['T']};
    qr[${nvars}] = 1.0;
%  for n,spn in enumerate(c['names'][:-1]):
    qr[${Yix+n}] = ${c[spn]};
    qr[${nvars}] -= ${c[spn]};
%  endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

    // We now have a valid density value, set momentums to reach input mdot
    fpdtype_t invrho = 1.0/ur[0];
% for i in range(ndims):
    ur[${i + 1}] = -2*nl[${i}]*${c['mdot-per-area']} - ul[${i+1}];
    qr[${i + 1}] = ur[${i+1}]*invrho;
% endfor

   // We have added ke, so we need to update total energy
    ur[${ndims+1}] += 0.5*invrho*${pyfr.dot('ur[{i}]', i=(1,ndims+1))};

</%pyfr:macro>