<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>
<% ns = c['ns'] %>
<% Yix = ndims+2 %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
    fpdtype_t invrho = 1.0/ul[0];

    // Set right primatives (except pressure)

    // Set zero velocity so we only compute internal energy
% for i in range(ndims):
    qr[${i+1}] = 0.0;
% endfor
    // Set Temperature
    qr[${ndims+1}] = ${c['T']};

    // Set species
    qr[${nvars}] = 1.0;
% for n in range(ns-1):
    qr[${Yix + n}] = ul[${Yix + n}]*rhoinv;
    qr[${nvars}] -= qr[${Yix + n}];
% endfor

    // Set right density
    ur[0] = ul[0];
    ${pyfr.expand('stateFrom-rhoTY', 'ur', 'qr', 'qhr')};

    // Have internal energy set, add momentum
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = -ul[${i + 1}] + 2*${c[v]}*ul[0];
    qr[${i+1}] = ur[${i+1}]*invrho;
% endfor

    // Add in the ke
    ur[${ndims + 1}] += 0.5*invrho*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
    fpdtype_t invrho = 1.0/ul[0];

    // Set right primatives (except pressure)

    // Set zero velocity so we only compute internal energy
% for i in range(ndims):
    qr[${i+1}] = 0.0;
% endfor
    // Set Temperature
    qr[${ndims+1}] = ${c['T']};

    // Set species
    qr[${nvars}] = 1.0;
% for n in range(ns-1):
    qr[${Yix + n}] = ul[${Yix + n}]*rhoinv;
    qr[${nvars}] -= qr[${Yix + n}];
% endfor

    // Set right density
    ur[0] = ul[0];
    ${pyfr.expand('stateFrom-rhoTY', 'ur', 'qr', 'qhr')};

    // Have internal energy set, add momentum
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = -ul[${i + 1}] + 2*${c[v]}*ul[0];
    qr[${i + 1}] = ur[${i + 1}]*invrho;
% endfor

    // Add in the ke
    ur[${ndims + 1}] += 0.5*invrho*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};

</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
