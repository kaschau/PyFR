<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = 'fast_coeff' in c %>\
<% niter_max = 4 if fast_props else 7 %>\
% if fast_props:
<% fast_coeff = c['fast_coeff'] %>\
% else:
<% T_cutoff = c['T_cutoff'] %>\
<% NASA7_Thigh = c['NASA7_Thigh'] %>\
<% NASA7_Tlow = c['NASA7_Tlow'] %>\
% endif\

<%pyfr:macro name='T_iter' params='e, cp, Rmix, T, q, qh'>

% for i in range(niter_max):
{
  fpdtype_t h = 0.0;
  cp = 0.0;
  fpdtype_t cpp = 0.0;
  % for n in range(ns):
  // ${c['names'][n]} Properties
  {
    fpdtype_t cps, hs, cpps;
    % if fast_props:
        cps = ${pyfr.nasa_cps(fast_coeff[n], Ru, MW[n])};
        hs = ${pyfr.nasa_hs(fast_coeff[n], Ru, MW[n])};
        cpps = ${pyfr.nasa_cpp(fast_coeff[n], Ru, MW[n])};
        cpp += cpps * q[${n}];
    % else:
      if (T < ${T_cutoff[n]})
      {
        cps = ${pyfr.nasa_cps(NASA7_Tlow[n], Ru, MW[n])};
        hs = ${pyfr.nasa_hs(NASA7_Tlow[n], Ru, MW[n])};
      }else
      {
        cps = ${pyfr.nasa_cps(NASA7_Thigh[n], Ru, MW[n])};
        hs = ${pyfr.nasa_hs(NASA7_Thigh[n], Ru, MW[n])};
      }
    % endif
    cp += cps * q[${n}];
    h += hs * q[${n}];
    % if i == niter_max - 1:
    qh[${4 + n}] = hs;
    % endif
  }
  % endfor
  fpdtype_t f = e - (h - Rmix * T);
  fpdtype_t fp = -cp + Rmix;
  % if fast_props:
  fpdtype_t fpp = -cpp;
  // Halleys's Method
  T -= (f*fp) / (fp*fp - 0.5*f*fpp);
  % else:
  // Newtons's Method
  T -= f / fp;
  % endif
}
% endfor
</%pyfr:macro>