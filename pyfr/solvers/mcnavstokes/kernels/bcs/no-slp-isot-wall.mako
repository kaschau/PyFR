<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    // Set right primatives
    qr[${pix}] = ql[${pix}];
    // Set wall velocity
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + vix}] = -ql[${i + vix}] + 2*${c[v]};
% endfor
    // Set Temperature
    qr[${Tix}] = ${c['T']};

    // Set species
% for n in range(ns):
    qr[${n}] = ql[${n}];
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    // Set right primatives
    qr[${pix}] = ql[${pix}];
    // Set wall velocity
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + vix}] = -ql[${i + vix}] + 2*${c[v]};
% endfor
    // Set Temperature
    qr[${Tix}] = ${c['T']};

    // Set species
% for n in range(ns):
    qr[${n}] = ql[${n}];
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, ql, qhl, nl, grad_ul, grad_ur'>
    fpdtype_t invrho = 1.0/ql[${rhoix}];

    // Copy non species fluid-side gradients across to wall-side gradients
% for i, j in pyfr.ndrange(ndims, ndims + 1):
    gradur[${i}][${j + vix}] = gradul[${i}][${j + vix}];
% endfor

% if ndims == 2:

    // Enforce zero normal species gradient in wall
    fpdtype_t rho_x = ${" + ".join([f"grad_ul[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_ul[1][{n}]" for n in range(ns)])};

    fpdtype_t Y;
%   for n in range(ns-1):
    Y = ql[${n}];
    grad_ur[0][${n}] = Y*rhol_x;
    grad_ur[1][${n}] = Y*rhol_y;
%   endfor

% elif ndims == 3:

    // Enforce zero normal species gradient in wall
    fpdtype_t rho_x = ${" + ".join([f"grad_ul[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_ul[1][{n}]" for n in range(ns)])};
    fpdtype_t rho_z = ${" + ".join([f"grad_ul[2][{n}]" for n in range(ns)])};

    fpdtype_t Y;
%   for n in range(ns-1):
    Y = ql[${n}];
    grad_ur[0][${n}] = Y*rhol_x;
    grad_ur[1][${n}] = Y*rhol_y;
    grad_ur[2][${n}] = Y*rhol_z;
%   endfor

% endif
</%pyfr:macro>