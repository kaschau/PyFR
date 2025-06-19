<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = 'fast_coeff' in c %>\
% if fast_props:
<% fast_coeff = c['fast_coeff'] %>\
% else:
<% T_cutoff = c['T_cutoff'] %>\
<% NASA7_Thigh = c['NASA7_Thigh'] %>\
<% NASA7_Tlow = c['NASA7_Tlow'] %>\
% endif\

<%pyfr:macro name='stateFrom-prims' params='u, q, qh'>

    ## q is an array of length nvars + 2
    ## storing all primitive s
    ## 0:ns-1,    ns:ns+ndims, ns+ndims+1, ns+ndims+2, nvars + 2
    ## Y0...Ynsp, u,v(,w),  rho          , p,          T

    ## qh stores mixture thermodynamic properties
    ## 0,  1,     2, 3, 4..4 + ns
    ## cp, gamma, c, e, hi1..hins


    // Compute mixture properties
    fpdtype_t R = 0.0;
% for n in range(ns):
    R += q[${n}]*${Ru/MW[n]};
% endfor

    // Compute density
    fpdtype_t rho = q[${pix}]/(R*q[${Tix}]);
    q[${rhoix}] = rho;

    // Species mass
% for n in range(ns):
    u[${n}] = q[${n}]*rho;
% endfor

    // Compute momentum
% for i in range(ndims):
    u[${i + vix}] = q[${i + vix}]*rho;
% endfor

    // Total energy
    fpdtype_t T = q[${Tix}];
    fpdtype_t h = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    // ${c['names'][n]} Properties
    {
      fpdtype_t cps, hs;
      % if fast_props:
          cps = ${pyfr.nasa_cps(fast_coeff[n], Ru, MW[n])};
          hs = ${pyfr.nasa_hs(fast_coeff[n], Ru, MW[n])};
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
      h += hs * q[${n}];
      cp += cps * q[${n}];
      qh[${4 + n}] = hs;
    }
% endfor

    fpdtype_t rhoe = rho*h - q[${pix}];
    u[${Eix}] = rhoe + 0.5*rho*${pyfr.dot('q[{i}]', i=(vix,vix + ndims))};

    // Store gamma, cp, c, rhoe
    qh[0] = cp / (cp - R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0]*R*q[${Tix}]);
    qh[3] = rhoe;

    // Store species enthalpy (per mass)
    // ^ done up there

</%pyfr:macro>