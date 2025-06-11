<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% MW = c['MW'] %>\
<% muPoly = c['muPoly'] %>\
<% kappaPoly = c['kappaPoly'] %>\
## Dij is a symmetric matrix, so just store half and index appropriately
<% DijPoly = c['DijPoly'] %>\
<% Dijix = lambda n,n2 : int(ns * (ns-1) / 2 - (ns - n) * (ns - n - 1)/2 + n2) %>\

<%pyfr:macro name='mixture_transport' params='u, q, qh, qt'>
  fpdtype_t invp = 1.0 / q[${pix}];
  fpdtype_t T = q[${Tix}];

  // Mole fraction
  fpdtype_t MWmix = 0.0;
  fpdtype_t X[${ns}];
  {
    fpdtype_t mass = 0.0;
% for n in range(ns):
    X[${n}] = q[${n}] * ${1.0 / MW[n]};
    mass += X[${n}];
% endfor
    fpdtype_t invmass = 1.0 / mass;
% for n in range(ns):
    X[${n}] *= invmass;
    MWmix += X[${n}] * ${MW[n]};
    // Avoid pure species condition
    X[${n}] = fmax(X[${n}], ${fpdtype_eps});
% endfor
  }

  // Evaluate all property poly'l
  fpdtype_t logT = log(T);
  fpdtype_t sqrtT = sqrt(T);
  fpdtype_t sqrtsqrtT = sqrt(sqrtT);
  fpdtype_t T_m3o2 = 1.0/(sqrtT*sqrtT*sqrtT);

  // Viscosity
  // because we fit poly'l to sqrt(mu), and we need sqrt(mu) a lot, just
  // evaluate the poly'l first, then convert to physical viscosity at the end
  fpdtype_t polymu_sp[${ns}];
  fpdtype_t invpolymu_sp[${ns}];
  % for n in range(ns):
    // ${c['names'][n]} viscosity
    <% deg = len(muPoly[n]) - 1%>\
    polymu_sp[${n}] = ${'+ logT*('.join(str(c) for c in muPoly[n])+')'*deg};
    invpolymu_sp[${n}] = 1.0/polymu_sp[${n}];
  % endfor
  {
    fpdtype_t mu = 0.0;
    fpdtype_t phitemp[${ns}] = {0};
    % for n in range(ns):
      // ${c['names'][n]} viscosity
      {
        % for n2 in range(n, ns):
          {
          % if n == n2:
            phitemp[${n}] += X[${n2}];
          % else:
            {
              fpdtype_t polymu_ratio = polymu_sp[${n}]*invpolymu_sp[${n2}];
              fpdtype_t num = 1.0 + polymu_ratio * ${(MW[n2] / MW[n])**0.25};
              fpdtype_t phi_temp = num * num * ${1.0/(math.sqrt(8.0) * math.sqrt(1.0 + MW[n]/MW[n2]))};
              phitemp[${n}] += phi_temp * X[${n2}];

              polymu_ratio = polymu_sp[${n2}]*invpolymu_sp[${n}];
              num = 1.0 + polymu_ratio * ${(MW[n] / MW[n2])**0.25};
              phi_temp = num * num * ${1.0/(math.sqrt(8.0)*math.sqrt(1.0 + MW[n2]/MW[n]))};
              phitemp[${n2}] += phi_temp * X[${n}];
            }
          % endif
          }
        % endfor
        // Set to correct dimensions
        fpdtype_t mu_sp = polymu_sp[${n}] * polymu_sp[${n}]*sqrtT;
        mu += mu_sp * X[${n}] / phitemp[${n}];
      }
    % endfor
    qt[0] = mu;
  }

  // Thermal conductivity
  {
    fpdtype_t sum1 = 0.0;
    fpdtype_t sum2 = 0.0;
    % for n in range(ns):
    {
      // ${c['names'][n]} thermal conductivity
      <% deg = len(kappaPoly[n]) - 1 %>\
      fpdtype_t kappa_sp = ${'+ logT*('.join(str(c) for c in kappaPoly[n])+')'*deg};
      kappa_sp *= sqrtT;
      sum1 += X[${n}] * kappa_sp;
      sum2 += X[${n}] / kappa_sp;
    }
    % endfor
    fpdtype_t kappa = 0.5*(sum1 + 1.0 / sum2);
    qt[1] = kappa;
  }

  // Diffusion coefficient
  {
    fpdtype_t sum1[${ns}] = {0};
    fpdtype_t sum2[${ns}] = {0};

    // Iterate over upper triangle only
    % for n in range(ns):
      % for n2 in range(n+1, ns):
      {
        <% ix = Dijix(n,n2) %>\
        <% deg = len(DijPoly[ix]) - 1 %>\
        fpdtype_t invDij = (${'+ logT*('.join(str(c) for c in DijPoly[ix])+')'*deg})*T_m3o2;
        // Contribute to species ${n} sums
        fpdtype_t temp = X[${n2}] * invDij;
        sum1[${n}] += temp;
        sum2[${n}] += temp * ${MW[n2]};

        // Contribute to species ${n2} sums (symmetric)
        temp = X[${n}] * invDij;
        sum1[${n2}] += temp;
        sum2[${n2}] += temp * ${MW[n]};
      }
      % endfor
    % endfor

    // Final computation for each species
    % for n in range(ns):
    {
      fpdtype_t final_sum2 = sum2[${n}] * X[${n}] / (MWmix - ${MW[n]} * X[${n}]);
      qt[${2 + n}] = invp / (sum1[${n}] + final_sum2);
    }
    % endfor
  }


#ifdef DEBUG
  printf("*********************************\n");
  printf("TRANSPORT PROPERTIES\n");
  printf("INPUT STATE\n");
  printf("trans&rho = %e\n", q[${rhoix}]);
  printf("trans&p = %e\n", q[${pix}]);
  printf("trans&T = %e\n", q[${Tix}]);
% for n in range(ns):
  printf("trans&Y_${c['names'][n]} = %e\n", q[${n}]);
% endfor

  printf("\nCOMPUTED PROPERTIES\n");
  printf("trans&MWmix = %e\n", MWmix);
  printf("trans&mu = %e\n", qt[0]);
  printf("trans&kappa = %e\n", qt[1]);
% for n in range(ns):
  printf("trans&D_${c['names'][n]} = %e\n", qt[${2+n}]);
% endfor
  printf("*********************************\n");

#endif
</%pyfr:macro>