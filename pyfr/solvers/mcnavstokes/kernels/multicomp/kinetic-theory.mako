<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%import math %>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='mix_visc_Wilke' params='polymu_sp, invpolymu_sp, X, sqrtT, mu'>
    fpdtype_t phitemp[${mcf.ns}] = {0};
    % for n in range(mcf.ns):
      {
        % for n2 in range(n, mcf.ns):
          {
          % if n == n2:
            phitemp[${n}] += X[${n2}];
          % else:
<%
            mwn = mcf[n].MW
            mwn2 = mcf[n2].MW
%>\
            {
              fpdtype_t polymu_ratio = polymu_sp[${n}]*invpolymu_sp[${n2}];
              fpdtype_t num = 1.0 + polymu_ratio * ${(mwn2 / mwn)**0.25};
              fpdtype_t phi_temp = num * num * ${1.0/(math.sqrt(8.0) * math.sqrt(1.0 + mwn/mwn2))};
              phitemp[${n}] += phi_temp * X[${n2}];

              polymu_ratio = polymu_sp[${n2}]*invpolymu_sp[${n}];
              num = 1.0 + polymu_ratio * ${(mwn / mwn2)**0.25};
              phi_temp = num * num * ${1.0/(math.sqrt(8.0)*math.sqrt(1.0 + mwn2/mwn))};
              phitemp[${n2}] += phi_temp * X[${n}];
            }
          % endif
          }
        % endfor
        fpdtype_t mu_sp = polymu_sp[${n}] * polymu_sp[${n}]*sqrtT;
        mu += mu_sp * X[${n}] / phitemp[${n}];
      }
    % endfor
</%pyfr:macro>

<%pyfr:macro name='mix_visc_Herning-Zipperer' params='polymu_sp, invpolymu_sp, X, sqrtT, mu'>
    fpdtype_t numerator = 0.0;
    fpdtype_t denominator = 0.0;
    % for n in range(mcf.ns):
    {
      fpdtype_t X_sqrt_MW = X[${n}] * ${math.sqrt(mcf[n].MW)};
      numerator += X_sqrt_MW * polymu_sp[${n}] * polymu_sp[${n}];
      denominator += X_sqrt_MW;
    }
    % endfor
    mu = (numerator / denominator) * sqrtT;
</%pyfr:macro>

<%pyfr:macro name='mixture_transport' params='u, q, qh, qt'>
  fpdtype_t p = q[${pix}];
  fpdtype_t T = q[${Tix}];

  // Mole fraction
  fpdtype_t MWmix = 0.0;
  fpdtype_t X[${mcf.ns}];
  {
    fpdtype_t total_moles = 0.0;
% for n in range(mcf.ns):
    X[${n}] = q[${n}] * ${1.0 / mcf[n].MW};
    total_moles += X[${n}];
% endfor
    fpdtype_t inv_total_moles = 1.0 / total_moles;
% for n in range(mcf.ns):
    X[${n}] *= inv_total_moles;
    MWmix += X[${n}] * ${mcf[n].MW};
    X[${n}] = fmax(X[${n}], ${fpdtype_eps});
% endfor
  }

  // Viscosity
  fpdtype_t logT = log(T);
  fpdtype_t sqrtT = sqrt(T);
  fpdtype_t polymu_sp[${mcf.ns}];
  fpdtype_t invpolymu_sp[${mcf.ns}];
  % for n in range(mcf.ns):
    polymu_sp[${n}] = ${mcf[n].mu_expr('logT')};
    invpolymu_sp[${n}] = 1.0/polymu_sp[${n}];
  % endfor

  {
  fpdtype_t mu = 0.0;
  ${pyfr.expand('mix_visc_' + mcf.mixing_rule, 'polymu_sp', 'invpolymu_sp', 'X', 'sqrtT', 'mu')};
  qt[0] = mu;
  }

  // Thermal conductivity
  {
    fpdtype_t sum1 = 0.0;
    fpdtype_t sum2 = 0.0;
    % for n in range(mcf.ns):
    {
      fpdtype_t kappa_sp = ${mcf[n].kappa_expr('logT')};
      kappa_sp *= sqrtT;
      sum1 += X[${n}] * kappa_sp;
      sum2 += X[${n}] / kappa_sp;
    }
    % endfor
    qt[1] = 0.5*(sum1 + 1.0 / sum2);
  }

  // Diffusion coefficient
  {
    fpdtype_t T_m3o2 = 1.0/(sqrtT*sqrtT*sqrtT);
    fpdtype_t sum1[${mcf.ns}] = {0};
    fpdtype_t sum2[${mcf.ns}] = {0};

    % for n in range(mcf.ns):
      % for n2 in range(n+1, mcf.ns):
      {
        fpdtype_t invDij = ${mcf[n].dij_expr(n2, 'logT')}*T_m3o2;
        fpdtype_t temp = X[${n2}] * invDij;
        sum1[${n}] += temp;
        sum2[${n}] += temp * ${mcf[n2].MW};

        temp = X[${n}] * invDij;
        sum1[${n2}] += temp;
        sum2[${n2}] += temp * ${mcf[n].MW};
      }
      % endfor
      sum2[${n}] *= X[${n}] / (MWmix - ${mcf[n].MW} * X[${n}]);
      qt[${2 + n}] = 1.0 / (p * (sum1[${n}] + sum2[${n}]));
    % endfor
  }


#ifdef DEBUG
  printf("*********************************\n");
  printf("TRANSPORT PROPERTIES\n");
  printf("INPUT STATE\n");
  printf("trans&rho = %e\n", q[${rhoix}]);
  printf("trans&p = %e\n", q[${pix}]);
  printf("trans&T = %e\n", q[${Tix}]);
% for n in range(mcf.ns):
  printf("trans&Y_${mcf.sp_names[n]} = %e\n", q[${n}]);
% endfor

  printf("\nCOMPUTED PROPERTIES\n");
  printf("trans&MWmix = %e\n", MWmix);
  printf("trans&mu = %e\n", qt[0]);
  printf("trans&kappa = %e\n", qt[1]);
% for n in range(mcf.ns):
  printf("trans&D_${mcf.sp_names[n]} = %e\n", qt[${2+n}]);
% endfor
  printf("*********************************\n");

#endif
</%pyfr:macro>
