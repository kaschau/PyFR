<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor

% for i in range(ndims):
    ur[${i + vix}] = -ul[${i + vix}];
% endfor
    ur[${Eix}] = ul[${Eix}];

    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor

% for i in range(ndims):
    ur[${i + vix}] = 0.0;
% endfor

    fpdtype_t rho = ${" + ".join([f"ul[{n}]" for n in range(ns)])};
    ur[${Eix}] = ul[${Eix}]
                     - (0.5/rho)*${pyfr.dot('ul[{i}]', i=(vix,vix + ndims))};
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, ql, qhl, nl, grad_ul, grad_ur'>
    fpdtype_t rhoE = ul[${Eix}];
    fpdtype_t rho = ql[${rhoix}];
    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t E = rhoE*rcprho;

    // Copy non species/energy fluid-side gradients across to wall-side gradients
% for i, j in pyfr.ndrange(ndims, ndims):
    gradur[${i}][${j+vix}] = gradul[${i}][${j+vix}];
% endfor

% if ndims == 2:

    fpdtype_t rho_x = ${" + ".join([f"grad_ul[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_ul[1][{n}]" for n in range(ns)])};

    // Velocity
    fpdtype_t u = ql[${vix}];
    fpdtype_t v = ql[${vix + 1}];

    // Velocity derivatives (rho*d[u,v]/d[x,y])
    fpdtype_t u_x = grad_ul[0][${vix}] - u*rho_x;
    fpdtype_t u_y = grad_ul[1][${vix}] - u*rho_y;
    fpdtype_t v_x = grad_ul[0][${vix + 1}] - v*rho_x;
    fpdtype_t v_y = grad_ul[1][${vix + 1}] - v*rho_y;

    // Enforce zero normal temperature gradient in wall
    fpdtype_t e_Y_Y_x = 0.0;
    fpdtype_t e_Y_Y_y = 0.0;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_y', 'ul', 'ql', 'qhl', 'grad_ul')};
    grad_ur[0][${Eix}] = u*u_x + v*v_x + rho*e_Y_Y_x + E*rho_x;
    grad_ur[1][${Eix}] = u*u_y + v*v_y + rho*e_Y_Y_y + E*rho_y;
    // Enforce zero normal species gradient in wall
%   for n in range(ns):
    grad_ur[0][${n}] = ql[${n}]*rho_x;
    grad_ur[1][${n}] = ql[${n}]*rho_y;
%   endfor

% elif ndims == 3:
    fpdtype_t rho_x = ${" + ".join([f"grad_ul[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_ul[1][{n}]" for n in range(ns)])};
    fpdtype_t rho_z = ${" + ".join([f"grad_ul[2][{n}]" for n in range(ns)])};

    // Velocity
    fpdtype_t u = ql[${vix}];
    fpdtype_t v = ql[${vix + 1}];
    fpdtype_t w = ql[${vix + 2}];

    // Velocity derivatives (rho*d[u,v]/d[x,y])
    fpdtype_t u_x = grad_ul[0][${vix}] - u*rho_x;
    fpdtype_t u_y = grad_ul[1][${vix}] - u*rho_y;
    fpdtype_t u_z = grad_ul[2][${vix}] - u*rho_z;
    fpdtype_t v_x = grad_ul[0][${vix + 1}] - v*rho_x;
    fpdtype_t v_y = grad_ul[1][${vix + 1}] - v*rho_y;
    fpdtype_t v_z = grad_ul[2][${vix + 1}] - v*rho_z;
    fpdtype_t w_x = grad_ul[0][${vix + 2}] - w*rho_x;
    fpdtype_t w_y = grad_ul[1][${vix + 2}] - w*rho_y;
    fpdtype_t w_z = grad_ul[2][${vix + 2}] - w*rho_z;

    // Enforce zero normal temperature gradient in wall
    fpdtype_t e_Y_Y_x = 0.0;
    fpdtype_t e_Y_Y_y = 0.0;
    fpdtype_t e_Y_Y_z = 0.0;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_y', 'e_Y_Y_z', 'ul', 'ql', 'qhl', 'grad_ul')};
    grad_ur[0][${Eix}] = u*u_x + v*v_x + rho*e_Y_Y_x + E*rho_x;
    grad_ur[1][${Eix}] = u*u_y + v*v_y + rho*e_Y_Y_y + E*rho_y;
    grad_ur[2][${Eix}] = u*u_z + v*v_z + rho*e_Y_Y_z + E*rho_z;
    // Enforce zero normal species gradient in wall
%   for n in range(ns):
    grad_ur[0][${n}] = ql[${n}]*rho_x;
    grad_ur[1][${n}] = ql[${n}]*rho_y;
    grad_ur[2][${n}] = ql[${n}]*rho_z;
%   endfor

% endif
</%pyfr:macro>
