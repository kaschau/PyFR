<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = -ul[${i + 1}];
% endfor
    ur[${ndims + 1}] = ul[${ndims + 1}];
<% ns = c['ns'] %>
% for n in range(ns-1):
    ur[${ndims + 2 + n}] = ul[${ndims + 2 + n}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = 0.0;
% endfor
    ur[${ndims+1}] = ul[${ndims+1}]
                     - (0.5/ul[0])*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))};
<% ns = c['ns'] %>
% for n in range(ns-1):
    ur[${ndims + 2 + n}] = ul[${ndims + 2 + n}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, nl, qhl, grad_ul, grad_ur'>
    fpdtype_t rhoE = ul[${ndims+1}];
    fpdtype_t rho = ul[0];
    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t E = rhoE*rcprho;

% if ndims == 2:
    fpdtype_t rhol_x = grad_ul[0][0];
    fpdtype_t rhol_y = grad_ul[1][0];

    // Copy all fluid-side gradients across to wall-side gradients
    ${pyfr.expand('bc_common_grad_copy', 'ul', 'nl', 'grad_ul', 'grad_ur')};

    // Enforce zero normal temperature gradient in wall
    grad_ur[0][3] = E*rhol_x;
    grad_ur[1][3] = E*rhol_y;
    // Enforce zero normal species gradient in wall
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>
    fpdtype_t Y;
%   for n in range(ns-1):
    Y = ul[${Yix+n}]*rcprho;
    grad_ur[0][${Yix+n}] = Y*rhol_x;
    grad_ur[1][${Yix+n}] = Y*rhol_y;
%   endfor

% elif ndims == 3:
    fpdtype_t rhol_x = grad_ul[0][0];
    fpdtype_t rhol_y = grad_ul[1][0];
    fpdtype_t rhol_z = grad_ul[2][0];

    // Copy all fluid-side gradients across to wall-side gradients
    ${pyfr.expand('bc_common_grad_copy', 'ul', 'nl', 'grad_ul', 'grad_ur')};

    // Enforce zero normal temperature gradient in wall
    grad_ur[0][3] = E*rhol_x;
    grad_ur[1][3] = E*rhol_y;
    grad_ur[2][3] = E*rhol_z;
    // Enforce zero normal species gradient in wall
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>
    fpdtype_t Y;
%   for n in range(ns-1):
    Y = ul[${Yix+n}]*rcprho;
    grad_ur[0][${Yix+n}] = Y*rhol_x;
    grad_ur[1][${Yix+n}] = Y*rhol_y;
    grad_ur[2][${Yix+n}] = Y*rhol_z;
%   endfor

% endif
</%pyfr:macro>
