<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.chem.net-rate-of-production'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\
<% MW = c['MW'] %>\
<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% nu_f = c['nu_f'] %>\
<% nu_b = c['nu_b'] %>\
<% fast_props = N7.shape[1] == 7 %>\
<% reconstruct = nsub_steps > 1 %>\
<% tSub = dt / float(nsub_steps) %>\

<%pyfr:macro name='finite_rate_source' params='t, u, ploc, src'>

  // Compute thermodynamic properties
  fpdtype_t q[${nvars + 2}];
  fpdtype_t qh[${4 + ns}];
  ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};

  fpdtype_t rho = q[${rhoix}];
  fpdtype_t rhoinv = 1.0/rho;
  fpdtype_t T = q[${Tix}];

  ## Generate straightforward finite rate source terms
  %if not reconstruct:

    ${pyfr.expand('net_rate_of_production', 'q', 'T', 'rho', 'src')};

  % else: ## take sub steps

  double Ynew[${ns}];
  double Yold[${ns}];
  double rhoY[${ns}];
  double tmpSrc[${ns}];
  % for n in range(ns):
    Yold[${n}] = q[${n}];
    rhoY[${n}] = u[${n}];
  % endfor
  for(int nSub = 0; nSub < ${nsub_steps}; nSub++){
    ${pyfr.expand('net_rate_of_production', 'Yold', 'T', 'rho', 'tmpSrc')};

    // Compute cp
    fpdtype_t cp = 0.0;
    % for n in range(ns):
    {
      % if fast_props:
      {
        fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 0)};
        cp += cps*Yold[${n}];
      }
      % else:
      if (T < ${N7[n,0]}){
        fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 8)};
        cp += cps*Yold[${n}];
      }else{
        fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 1)};
        cp += cps*Yold[${n}];
      }
      % endif
    }
    % endfor

    // Take sub step in time for species
    fpdtype_t Yact_sum = 0.0;
    fpdtype_t Ybath_sum = 0.0;
    % for n in range(ns):
    {
      <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
        Ynew[${n}] = Yold[${n}] + tmpSrc[${n}] * rhoinv * ${tSub};
        Ynew[${n}] = fmax(0.0, Ynew[${n}]);
        Yact_sum += Ynew[${n}];
      % else:
        Ynew[${n}] = Yold[${n}];
        Ybath_sum += Ynew[${n}];
      % endif
    }
    % endfor
    // Normalize the active species (non-bath) and their sources
    fpdtype_t Y_norminv = (1.0-Ybath_sum)/Yact_sum;
    % for n in range(ns):
      <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
      Ynew[${n}] *= Y_norminv;
      tmpSrc[${n}] = rho*(Ynew[${n}] - Yold[${n}]) * ${1.0/dt};
      % endif
    % endfor

    // Take sub step in time for temperature
    fpdtype_t dTdt = 0.0;
    fpdtype_t Tinv = 1.0/T;
    % for n in range(ns):
    {
      <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
        % if fast_props:
          fpdtype_t hi = ${pyfr.nasa_hi(N7[n,:], 0)};
        % else:
          fpdtype_t hi;
          if (T < ${N7[n,0]})
          {
            hi = ${pyfr.nasa_hi(N7[n,:], 8)};
          }else
          {
            hi = ${pyfr.nasa_hi(N7[n,:], 1)};
          }
        % endif
        dTdt -= hi * tmpSrc[${n}];
      % endif
    Yold[${n}] = Ynew[${n}];
    }
    % endfor
    dTdt /= cp * rho;
    T += dTdt * ${tSub};
  }

  // Reconstruct d(rhoY)/dt based on where we ended up
  % for n in range(ns):
    // ${c['names'][n]}
    <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
    % if max(abs(nu_sum)) > 0.0:
        src[${n}] = (Ynew[${n}] * rho - rhoY[${n}]) * ${1.0/dt};
    % else:
      src[${n}] = 0.0;
    % endif
  % endfor

% endif


// Set non chemical terms to zero
% for i in range(ndims):
  src[${i + vix}] = 0.0;
% endfor
  src[${Eix}] = 0.0;


#ifdef DEBUG
  printf("*********************************\n");
  printf("CHEMICAL SOURCE TERMS\n");
% for n in range(ns):
  printf("chem&omega_${c['names'][n]} = %e\n", src[${n}]);
% endfor
  printf("*********************************\n");
#endif

</%pyfr:macro>