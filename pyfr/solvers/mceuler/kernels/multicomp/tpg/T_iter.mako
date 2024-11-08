<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = N7.shape[1] == 7 %>\
<% niter_max = 4 if fast_props else 6 %>\

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
        cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 0)};
        hs = ${pyfr.nasa_hs(N7[n,:], Ru, MW[n], 0)};
        cpps = ${pyfr.nasa_cpp(N7[n,:], Ru, MW[n], 0)};
        cpp += cpps * q[${n}];
    % else:
      if (T < ${N7[n,0]})
      {
        cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 8)};
        hs = ${pyfr.nasa_hs(N7[n,:], Ru, MW[n], 8)};
      }else
      {
        cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 1)};
        hs = ${pyfr.nasa_hs(N7[n,:], Ru, MW[n], 1)};
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