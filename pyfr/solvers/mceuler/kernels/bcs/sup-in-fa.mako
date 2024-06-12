<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.prims_to_cons'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>

    fpdtype_t q[{nvars+1}];
    q[0] = ${c['p']};
% for i, v in enumerate('uvw'[:ndims]):
    q[${i + 1}] = ${c[v]};
% endfor
    q[${ndims+1}] = ${c['T']};
<% Yix = ndims + 2 %>
    q[${nvars}] = 1.0
%  for n,spn in enumerate(props['names'][:-1]):
    q[${Yix+n}] = ${c[spn]};
    q[${nvars}] -= ${c[spn]};
%  endfor

    ${pyfr.expand('prims_to_cons', 'q', 'u')};

</%pyfr:macro>
