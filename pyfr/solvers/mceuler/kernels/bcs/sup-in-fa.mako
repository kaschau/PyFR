<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.prims_to_cons'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>

    fpdtype_t qr[${nvars+1}];
    qr[0] = ${c['p']};
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + 1}] = ${c[v]};
% endfor
    qr[${ndims+1}] = ${c['T']};
<% Yix = ndims + 2 %>
    qr[${nvars}] = 1.0;
%  for n,spn in enumerate(c['names'][:-1]):
    qr[${Yix+n}] = ${c[spn]};
    qr[${nvars}] -= ${c[spn]};
%  endfor

    ${pyfr.expand('prims_to_cons', 'qr', 'ur')};

</%pyfr:macro>
