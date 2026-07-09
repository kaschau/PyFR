<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ${fluid.decl('ul', 'rho, invrho, v, Y, T, p, R, gammamix, a', suffix='wc')}

    // Outer (freestream) state from the config
    fpdtype_t Yb[${ns}];
% for n, spn in enumerate(sp_names):
    Yb[${n}] = ${c[spn]};
% endfor
    fpdtype_t Rb = ${' + '.join(f'Yb[{n}]*{r}'
                                for n, r in enumerate(sp_R))};
    fpdtype_t cpb = ${' + '.join(f'Yb[{n}]*{cp}'
                                 for n, cp in enumerate(sp_cp0))};
    fpdtype_t gb_o = cpb/(cpb - Rb);
    fpdtype_t pb_o = ${c['p']};
    fpdtype_t Tb_o = ${c['T']};
    fpdtype_t rhob_o = pb_o/(Rb*Tb_o);
    fpdtype_t cb_o = sqrt(gb_o*Rb*Tb_o);

    // Entropy and Riemann invariant ratios for the two states
    fpdtype_t s_o = pb_o*pow(rhob_o, -gb_o);
    fpdtype_t ratio_o = 2.0*cb_o/(gb_o - 1.0);
    fpdtype_t s_i = pwc*pow(rhowc, -gammamixwc);
    fpdtype_t ratio_i = 2.0*awc/(gammamixwc - 1.0);

    // Normal velocities of the outer and interior states
    fpdtype_t V_e = ${' + '.join(f"({c['uvw'[i]]})*nl[{i}]"
                                 for i in range(ndims))};
    fpdtype_t V_i = ${' + '.join(f'vwc[{i}]*nl[{i}]'
                                 for i in range(ndims))};

    fpdtype_t R_e = (fabs(V_e) >= cb_o && V_i >= 0)
                  ? V_i - ratio_i
                  : V_e - ratio_o;
    fpdtype_t gmo_e = (fabs(V_e) >= cb_o && V_i >= 0)
                    ? gammamixwc - 1.0
                    : gb_o - 1.0;
    fpdtype_t R_i = (fabs(V_e) >= cb_o && V_i < 0)
                  ? V_e + ratio_o
                  : V_i + ratio_i;
    fpdtype_t gmo_i = (fabs(V_e) >= cb_o && V_i < 0)
                    ? gb_o - 1.0
                    : gammamixwc - 1.0;

    // Boundary values
    fpdtype_t V_b = 0.5*(R_e + R_i);
    fpdtype_t c_b = 0.25*(R_i*gmo_i - R_e*gmo_e);
    fpdtype_t g_b = (V_i < 0) ? gb_o : gammamixwc;
    fpdtype_t R_b = (V_i < 0) ? Rb : Rwc;
    fpdtype_t s_b = (V_i < 0) ? s_o : s_i;
    fpdtype_t T_b = c_b*c_b/(g_b*R_b);
    fpdtype_t p_b = pow(s_b/pow(R_b*T_b, g_b), 1.0/(1.0 - g_b));

    fpdtype_t vb[${ndims}];
% for i in range(ndims):
    vb[${i}] = (V_i >= 0)
             ? vwc[${i}] + (V_b - V_i)*nl[${i}]
             : (${c['uvw'[i]]}) + (V_b - V_e)*nl[${i}];
% endfor

    // Assemble the boundary state (freestream composition)
    fpdtype_t rhob = p_b/(Rb*T_b);
    fpdtype_t hb = ${' + '.join(f"Yb[{n}]*({fluid.h_of_T(n, 'T_b')})"
                                for n in range(ns))};

% for n in range(ns):
    ur[${n}] = rhob*Yb[${n}];
% endfor
% for i in range(ndims):
    ur[${ns + i}] = rhob*vb[${i}];
% endfor
    ur[${nvars - 1}] = rhob*hb - p_b
        + 0.5*rhob*(${' + '.join(f'vb[{i}]*vb[{i}]'
                                 for i in range(ndims))});
</%pyfr:macro>
