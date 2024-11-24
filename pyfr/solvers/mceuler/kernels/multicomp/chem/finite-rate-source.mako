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

  for(int nSub = 0; nSub < ${nsub_steps}; nSub++){
    ${pyfr.expand('net_rate_of_production', 'q', 'T', 'rho', 'src')};

    // Compute cp
    fpdtype_t cp = 0.0;
    % for n in range(ns):
    {
      % if fast_props:
      {
        fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 0)};
        cp += cps*q[${n}];
      }
      % else:
      if (T < ${N7[n,0]}){
        fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 8)};
        cp += cps*q[${n}];
      }else{
        fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 1)};
        cp += cps*q[${n}];
      }
      % endif
    }
    % endfor

    // Take sub step in time
    fpdtype_t dTdt = 0.0;
    fpdtype_t tempsum = 0.0;
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
        dTdt -= hi * src[${n}];
        q[${n}] += src[${n}] * rhoinv * ${tSub};
        q[${n}] = fmax(0.0, q[${n}]);
      % endif
      tempsum += q[${n}];
    }
    % endfor
    // Normalize
    fpdtype_t tempsuminv = 1.0/tempsum;
    % for n in range(ns):
      q[${n}] *= tempsuminv;
    % endfor
    dTdt /= cp * rho;
    T += dTdt * ${tSub};
  }

  // Reconstruct d(rhoY)/dt based on where we ended up
  % for n in range(ns):
    // ${c['names'][n]}
    <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
    % if max(abs(nu_sum)) > 0.0:
        src[${n}] = (q[${n}] * rho - u[${n}]) * ${1.0/dt};
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