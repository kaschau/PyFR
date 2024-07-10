<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% ns = c['ns'] %>
<% Yix = ndims+2 %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
    fpdtype_t invrho = 1.0/ul[0];

    // Set right primatives
    qr[0] = ql[0];
    // Set wall velocity
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + 1}] = -ql[${i + 1}] + 2*${c[v]};
% endfor
    // Set Temperature
    qr[${ndims+1}] = ${c['T']};

    // Set species
    qr[${nvars}] = 1.0;
% for n in range(ns-1):
    qr[${Yix + n}] = ul[${Yix + n}]*invrho;
    qr[${nvars}] -= qr[${Yix + n}];
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
    fpdtype_t invrho = 1.0/ul[0];

    // Set right primatives
    qr[0] = ql[0];
    // Set wall velocity
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + 1}] = -ql[${i + 1}] + 2*${c[v]};
% endfor
    // Set Temperature
    qr[${ndims+1}] = ${c['T']};

    // Set species
    qr[${nvars}] = 1.0;
% for n in range(ns-1):
    qr[${Yix + n}] = ul[${Yix + n}]*invrho;
    qr[${nvars}] -= qr[${Yix + n}];
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, ql, qhl, nl, grad_ul, grad_ur'>
    fpdtype_t rcprho = 1.0/ul[0];

    // Copy non species fluid-side gradients across to wall-side gradients
% for i, j in pyfr.ndrange(ndims, ndims+2):
    gradur[${i}][${j}] = gradul[${i}][${j}];
% endfor

% if ndims == 2:
    // Enforce zero normal species gradient in wall
    fpdtype_t rhol_x = grad_ul[0][0];
    fpdtype_t rhol_y = grad_ul[1][0];
    fpdtype_t Y;
%   for n in range(ns-1):
    Y = ul[${Yix+n}]*rcprho;
    grad_ur[0][${Yix+n}] = Y*rhol_x;
    grad_ur[1][${Yix+n}] = Y*rhol_y;
%   endfor
% elif ndims == 3:
    // Enforce zero normal species gradient in wall
    fpdtype_t Y;
    fpdtype_t rhol_x = grad_ul[0][0];
    fpdtype_t rhol_y = grad_ul[1][0];
    fpdtype_t rhol_z = grad_ul[2][0];
%   for n in range(ns-1):
    Y = ul[${Yix+n}]*rcprho;
    grad_ur[0][${Yix+n}] = Y*rhol_x;
    grad_ur[1][${Yix+n}] = Y*rhol_y;
    grad_ur[2][${Yix+n}] = Y*rhol_z;
%   endfor
% endif
</%pyfr:macro>