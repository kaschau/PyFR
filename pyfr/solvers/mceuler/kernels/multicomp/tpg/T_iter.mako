<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='T_iter_newton' params='e, cp, Rmix, T, q, qh'>
    T = ${0.5*(mcf[0].thermo_ranges[0] + mcf[0].thermo_ranges[-1])};
% for i in range(mcf.T_iter_count):
{
  fpdtype_t h = 0.0;
  cp = 0.0;
  % for n in range(mcf.ns):
  {
    fpdtype_t cps = ${mcf[n].cp_expr('T')};
    fpdtype_t hs = ${mcf[n].h_expr('T')};
    cp += cps * q[${n}];
    h += hs * q[${n}];
  }
  % endfor
  T -= (e - (h - Rmix * T)) / (-cp + Rmix);
}
% endfor
    // Evaluate the thermodynamic state at the final temperature
    cp = 0.0;
% for n in range(mcf.ns):
{
    fpdtype_t cps = ${mcf[n].cp_expr('T')};
    fpdtype_t hs = ${mcf[n].h_expr('T')};
    cp += cps * q[${n}];
    qh[${4 + n}] = hs;
}
% endfor
</%pyfr:macro>

<%pyfr:macro name='T_iter_halley' params='e, cp, Rmix, T, q, qh'>
<%
    tc = [(sp.thermo_coeffs[0], mcf.Ru / sp.MW) for sp in mcf.species]
%>\
    fpdtype_t _a = -(${'+'.join([f'{c[1]*s/2.0}*q[{n}]' for n, (c, s) in enumerate(tc)])});
    fpdtype_t _b = Rmix - (${'+'.join([f'{c[0]*s}*q[{n}]' for n, (c, s) in enumerate(tc)])});
    fpdtype_t _c = e - (${'+'.join([f'{c[-2]*s}*q[{n}]' for n, (c, s) in enumerate(tc)])});
    T = fmin(${mcf[0].thermo_ranges[1]}, fmax(${mcf[0].thermo_ranges[0]}, fabs(_a) < ${fpdtype_eps} ? -_c/_b : (-_b + sqrt(fmax(0.0,_b*_b-4*_a*_c)))/(2*_a)));
% for i in range(mcf.T_iter_count):
{
  fpdtype_t h = 0.0;
  cp = 0.0;
  fpdtype_t cpp = 0.0;
  % for n in range(mcf.ns):
  {
    fpdtype_t cps = ${mcf[n].cp_expr('T')};
    fpdtype_t hs = ${mcf[n].h_expr('T')};
    fpdtype_t cpps = ${mcf[n].dcp_expr('T')};
    cp += cps * q[${n}];
    h += hs * q[${n}];
    cpp += cpps * q[${n}];
  }
  % endfor
  fpdtype_t f = e - (h - Rmix * T);
  fpdtype_t fp = -cp + Rmix;
  T -= (f*fp) / (fp*fp - 0.5*f*(-cpp));
}
% endfor
    // Evaluate the thermodynamic state at the final temperature
    cp = 0.0;
% for n in range(mcf.ns):
{
    fpdtype_t cps = ${mcf[n].cp_expr('T')};
    fpdtype_t hs = ${mcf[n].h_expr('T')};
    cp += cps * q[${n}];
    qh[${4 + n}] = hs;
}
% endfor
</%pyfr:macro>
