<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    // set right side primatives
%  for n,spn in enumerate(c['names']):
    qr[${n}] = ${c[spn]};
%  endfor

% for i in range(ndims):
    qr[${i + vix}] = ${c['uvw'[i]]};
% endfor

    qr[${pix}] = ${c['p']};
    qr[${Tix}] = ${c['T']};

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

    fpdtype_t gammar = qhr[0];
    fpdtype_t gammal = qhl[0];

    fpdtype_t gmor = gammar - 1.0;
    fpdtype_t gmol = gammal - 1.0;

    fpdtype_t cs = qhr[2];
    fpdtype_t s = ${c['p']}*pow(qr[${rhoix}], -gammar);
    fpdtype_t ratio = cs*2.0/gmor;

    fpdtype_t inv = 1.0/ql[${rhoix}];
    fpdtype_t V_e = ${' + '.join('{0}*nl[{1}]'.format(c['uvw'[i]], i)
                                 for i in range(ndims))};

    fpdtype_t V_i = ${' + '.join(f'ql[{vix + i}]*nl[{i}]' for i in range(ndims))};
    fpdtype_t p_i = ql[${pix}];
    fpdtype_t c_i = qhl[2];
    fpdtype_t R_e = (fabs(V_e) >= cs && V_i >= 0)
                  ? V_i - c_i*2.0/gmol
                  : V_e - ratio;
    fpdtype_t R_i = (fabs(V_e) >= cs && V_i < 0)
                  ? V_e + ratio
                  : V_i + c_i*2.0/gmol;
    fpdtype_t V_b = 0.5*(R_e + R_i);
    fpdtype_t c_b = 0.25*0.5*(gmor+gmol)*(R_i - R_e);
    fpdtype_t rho_b = (V_i < 0)
                    ? pow((1.0/(gammal*s))*c_b*c_b, 1.0/gmol)
                    : ql[${rhoix}]*pow(ql[${rhoix}]*c_b*c_b/(gammal*p_i), 1.0/gmol);
    fpdtype_t p_b = 1.0/gammal*rho_b*c_b*c_b;

%  for n,spn in enumerate(c['names']):
    ur[${n}] = ${c[spn]}*rho_b;
%  endfor

% for i in range(ndims):
    ur[${vix + i}] = (V_i >= 0)
                 ? rho_b*(ul[${vix + i}]*inv + (V_b - V_i)*nl[${i}])
                 : rho_b*(${c['uvw'[i]]} + (V_b - V_e)*nl[${i}]);
% endfor
    ur[${Eix}] = p_b*1.0/gmor
                     + 0.5*(1.0/rho_b)*${pyfr.dot('ur[{i}]', i=(vix,vix + ndims))};

    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>
