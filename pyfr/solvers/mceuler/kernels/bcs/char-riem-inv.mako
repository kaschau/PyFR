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

    fpdtype_t gr = qhr[0];
    fpdtype_t gl = qhl[0];

    fpdtype_t gmor = gr - 1.0;
    fpdtype_t gmol = gl - 1.0;

    fpdtype_t c_o = qhr[2];
    fpdtype_t s_o = ${c['p']}*pow(qr[${rhoix}], -gr);
    fpdtype_t ratio_o = c_o*2.0/gmor;

    fpdtype_t c_i = qhl[2];
    fpdtype_t s_i = ql[${pix}]*pow(ql[${rhoix}], -gl);
    fpdtype_t ratio_i = c_i*2.0/gmol;

    fpdtype_t V_e = ${' + '.join('{0}*nl[{1}]'.format(c['uvw'[i]], i)
                                 for i in range(ndims))};

    fpdtype_t V_i = ${' + '.join(f'ql[{vix + i}]*nl[{i}]' for i in range(ndims))};
    fpdtype_t p_i = ql[${pix}];
    fpdtype_t R_e = (fabs(V_e) >= c_o && V_i >= 0)
                  ? V_i - ratio_i
                  : V_e - ratio_o;
    fpdtype_t gmo_e = (fabs(V_e) >= c_o && V_i >= 0)
                    ? gmol
                    : gmor;
    fpdtype_t R_i = (fabs(V_e) >= c_o && V_i < 0)
                  ? V_e + ratio_o
                  : V_i + ratio_i;
    fpdtype_t gmo_i = (fabs(V_e) >= c_o && V_i < 0)
                    ? gmor
                    : gmol;

    // boundary values
    fpdtype_t V_b = 0.5*(R_e + R_i);
    fpdtype_t c_b = 0.25*(R_i*gmo_i - R_e*gmo_e);
    fpdtype_t g_b = (V_i < 0)
                  ? gr
                  : gl;
    fpdtype_t R_b = (V_i < 0)
                  ? qhr[1]*(gmor/gr)
                  : qhl[1]*(gmol/gl);
    fpdtype_t s_b = (V_i < 0)
                  ? s_o
                  : s_i;
    fpdtype_t T_b = c_b*c_b/(g_b*R_b);
    fpdtype_t p_b = pow(s_b/pow(R_b*T_b, g_b), 1.0/(1.0-g_b));

    qr[${pix}] = p_b;
    qr[${Tix}] = T_b;

% for i in range(ndims):
    qr[${vix + i}] = (V_i >= 0)
                 ? (ql[${vix + i}] + (V_b - V_i)*nl[${i}])
                 : (${c['uvw'[i]]} + (V_b - V_e)*nl[${i}]);
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};
</%pyfr:macro>
