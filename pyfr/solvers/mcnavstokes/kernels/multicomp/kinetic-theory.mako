<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% MW = c['MW'] %>\
<% muPoly = c['muPoly'] %>\
<% kappaPoly = c['kappaPoly'] %>\
## Dij is a symmetric matrix, so just store half and index appropriately
<% DijPoly = c['DijPoly'] %>\
<% Dijix = lambda n,n2 : int(ns * (ns-1) / 2 - (ns - n) * (ns - n - 1)/2 + n2) %>\

<%pyfr:macro name='mixture_transport' params='u, q, qh, qt'>
  fpdtype_t p = q[${pix}];
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

  fpdtype_t mu_sp[${ns}];
  fpdtype_t mu_spinv[${ns}];
  fpdtype_t invDij[${int((ns + 1)*ns/2)}];
% for n in range(ns):
  // ${c['names'][n]} viscosity, diffusion coefficients
  <% deg = len(muPoly[n]) - 1%>\
  mu_sp[${n}] = ${'+ logT*('.join(str(c) for c in muPoly[n])+')'*deg};

  // Set to correct dimensions
  mu_sp[${n}] *= sqrtsqrtT;
  mu_sp[${n}] *= mu_sp[${n}];
  mu_spinv[${n}] = 1.0/mu_sp[${n}];

  // Dont need to store every kappa!!!

  % for n2 in range(n, ns):
    <% ix = Dijix(n,n2)%>\
    <% deg = len(DijPoly[ix]) - 1 %>\
    invDij[${ix}] = (${'+ logT*('.join(str(c) for c in DijPoly[ix])+')'*deg})*T_m3o2;
  % endfor
% endfor

  // Now every species' property is computed, generate mixture values
  fpdtype_t mu = 0.0;
  fpdtype_t sum1 = 0.0;
  fpdtype_t sum2 = 0.0;

% for n in range(ns):
  // ${c['names'][n]} viscosity
  {
    fpdtype_t phitemp = 0.0;
    fpdtype_t sumd1 = 0.0;
    fpdtype_t sumd2 = 0.0;
    % for n2 in range(ns):
      {
        fpdtype_t num = 1.0 + sqrt(mu_sp[${n}] * mu_spinv[${n2}] * ${math.sqrt(MW[n2] / MW[n])});
        fpdtype_t phi = num*num*${1.0/(math.sqrt(8.0) * math.sqrt(1.0 + MW[n]/MW[n2]))};
        phitemp += phi * X[${n2}];
      % if n != n2:
      ##Symmetric
      <% ix = Dijix(n,n2) if n2>=n else Dijix(n2,n)%>\
        sumd1 += X[${n2}] * invDij[${ix}];
        sumd2 += X[${n2}] * ${MW[n2]} * invDij[${ix}];
      % endif
      }
    % endfor
    {
    mu += mu_sp[${n}] * X[${n}] / phitemp;
    // Mixture species diffusion coefficient
    // account for pressure
    sumd1 *= p;
    sumd2 *= p * X[${n}] / (MWmix - ${MW[n]} * X[${n}]);
    qt[${2 + n}] = 1.0 / (sumd1 + sumd2);

    // ${c['names'][n]} thermal conductivity
    <% deg = len(kappaPoly[n]) - 1 %>\
    fpdtype_t kappa_sp = ${'+ logT*('.join(str(c) for c in kappaPoly[n])+')'*deg};
    kappa_sp *= sqrtT;
    sum1 += X[${n}] * kappa_sp;
    sum2 += X[${n}] / kappa_sp;
    }
  }
% endfor

fpdtype_t kappa = 0.5*(sum1 + 1.0 / sum2);
qt[1] = kappa;
qt[0] = mu;


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