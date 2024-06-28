<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.prims_to_cons'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
<% ns = c['ns'] %>
<% Yix = ndims+2 %>

    // interior mixture state
    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('mixture_state', 'ul', 'ql', 'qhl')};

    // set right side primatives
    fpdtype_t qr[${nvars+1}];
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

    ${pyfr.expand('prims_to_cons', 'qr', 'ur')};

    // We now have a valid density value, set momentums to reach
    //input mdot
% for i in range(ndims):
    ur[${i + 1}] = nl[${i}]*${c['mdot-per-area']};
% endfor
   // We have added ke, so we need to update total energy
    ur[${ndims+1}] += 0.5/ur[0]*${pyfr.dot('ur[{i}]', i=(1,ndims+1))};

</%pyfr:macro>