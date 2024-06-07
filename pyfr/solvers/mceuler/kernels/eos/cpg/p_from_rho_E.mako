<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='p_from_rho_E' params='p, E, invrho, rhov'>
    p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)});
</%pyfr:macro>
