<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns = c['ns'] %>
<% Yix = ndims+2 %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = -ul[${i + 1}];
% endfor
    ur[${ndims + 1}] = ul[${ndims + 1}];
% for n in range(ns-1):
    ur[${Yix + n}] = ul[${Yix + n}];
% endfor

    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = 0.0;
% endfor
    ur[${ndims+1}] = ul[${ndims+1}]
                     - (0.5/ul[0])*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))};
% for n in range(ns-1):
    ur[${Yix + n}] = ul[${Yix + n}];
% endfor
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, ql, qhl, nl, grad_ul, grad_ur'>
    fpdtype_t rhoE = ul[${ndims+1}];
    fpdtype_t rho = ul[0];
    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t E = rhoE*rcprho;

    // Copy non energy/species fluid-side gradients across to wall-side gradients
% for i, j in pyfr.ndrange(ndims, ndims+1):
    gradur[${i}][${j}] = gradul[${i}][${j}];
% endfor

% if ndims == 2:
    fpdtype_t rho_x = grad_ul[0][0];
    fpdtype_t rho_y = grad_ul[1][0];

    // Velocity
    fpdtype_t u = ul[1]*rcprho;
    fpdtype_t v = ul[2]*rcprho;

    // Velocity derivatives (rho*d[u,v]/d[x,y])
    fpdtype_t u_x = grad_ul[0][1] - u*rho_x;
    fpdtype_t u_y = grad_ul[1][1] - u*rho_y;
    fpdtype_t v_x = grad_ul[0][2] - v*rho_x;
    fpdtype_t v_y = grad_ul[1][2] - v*rho_y;

    // Enforce zero normal temperature gradient in wall
    fpdtype_t e_Y_Y_x = 0.0;
    fpdtype_t e_Y_Y_y = 0.0;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_y', 'ul', 'ql', 'qhl', 'grad_ul')};
    grad_ur[0][3] = u*u_x + v*v_x + rho*e_Y_Y_x + E*rho_x;
    grad_ur[1][3] = u*u_y + v*v_y + rho*e_Y_Y_y + E*rho_y;
    // Enforce zero normal species gradient in wall
    fpdtype_t Y;
%   for n in range(ns-1):
    Y = ul[${Yix+n}]*rcprho;
    grad_ur[0][${Yix+n}] = Y*rho_x;
    grad_ur[1][${Yix+n}] = Y*rho_y;
%   endfor

% elif ndims == 3:
    fpdtype_t rho_x = grad_ul[0][0];
    fpdtype_t rho_y = grad_ul[1][0];
    fpdtype_t rho_z = grad_ul[2][0];

    // Velocity
    fpdtype_t u = ul[1]*rcprho;
    fpdtype_t v = ul[2]*rcprho;
    fpdtype_t w = ul[3]*rcprho;

    // Velocity derivatives (rho*grad[u,v,w])
    fpdtype_t u_x = grad_ul[0][1] - u*rho_x;
    fpdtype_t u_y = grad_ul[1][1] - u*rho_y;
    fpdtype_t u_z = grad_ul[2][1] - u*rho_z;
    fpdtype_t v_x = grad_ul[0][2] - v*rho_x;
    fpdtype_t v_y = grad_ul[1][2] - v*rho_y;
    fpdtype_t v_z = grad_ul[2][2] - v*rho_z;
    fpdtype_t w_x = grad_ul[0][3] - w*rho_x;
    fpdtype_t w_y = grad_ul[1][3] - w*rho_y;
    fpdtype_t w_z = grad_ul[2][3] - w*rho_z;

    // Enforce zero normal temperature gradient in wall
    fpdtype_t e_Y_Y_x = 0.0;
    fpdtype_t e_Y_Y_y = 0.0;
    fpdtype_t e_Y_Y_z = 0.0;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_y', 'e_Y_Y_z', 'ul', 'ql', 'qhl', 'grad_ul')};
    grad_ur[0][3] = u*u_x + v*v_x + w*w_x + rho*e_Y_Y_x + E*rho_x;
    grad_ur[1][3] = u*u_y + v*v_y + w*w_y + rho*e_Y_Y_y + E*rho_y;
    grad_ur[2][3] = u*u_z + v*v_z + w*w_z + rho*e_Y_Y_z + E*rho_z;

    // Enforce zero normal species gradient in wall
    fpdtype_t Y;
%   for n in range(ns-1):
    Y = ul[${Yix+n}]*rcprho;
    grad_ur[0][${Yix+n}] = Y*rho_x;
    grad_ur[1][${Yix+n}] = Y*rho_y;
    grad_ur[2][${Yix+n}] = Y*rho_z;
%   endfor

% endif
</%pyfr:macro>
