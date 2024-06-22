<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='finite_rate_source' params='t, u, ploc, src'>
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>

  fpdtype_t rho = u[0];


  // Output source terms
  src[0] = 0.0;
% for i in range(ndims):
  src[${i}] = 0.0;
% endfor
  src[${ndims+1}] = 0.0;
% for n in range(ns-1):
  src[${Yix+n}] += omega[${n}];
% endfor

</%pyfr:macro>