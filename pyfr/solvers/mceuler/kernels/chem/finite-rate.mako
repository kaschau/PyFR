<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.chem.net-rate-of-production'/>

<%pyfr:macro name='finite_rate' params='t, u, ploc, src'>
    ${fluid.decl('u', 'rho, Y, T', suffix='ch')}

    ${pyfr.expand('net_rate_of_production', 'Ych', 'Tch', 'rhoch', 'src')};

% for i in range(ndims):
    src[${ns + i}] = 0.0;
% endfor
    src[${nvars - 1}] = 0.0;
</%pyfr:macro>
