<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='rsolve_1d' params='ul, ur, ql, qr, qhl, qhr, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr, fsl, fsr;
    fpdtype_t usl[${nvars}], usr[${nvars}];

    ${pyfr.expand('inviscid_flux_1d', 'ul', 'fl', 'ql')};
    ${pyfr.expand('inviscid_flux_1d', 'ur', 'fr', 'qr')};

    % for i in range(ndims):
      vl[${i}] = ql[${i + vix}];
      vr[${i}] = qr[${i + vix}];
    % endfor
    pl = ql[${pix}];
    pr = qr[${pix}];

    fpdtype_t sqrtrl = sqrt(ql[${rhoix}]);
    fpdtype_t sqrtrr = sqrt(qr[${rhoix}]);
    fpdtype_t rr = qr[${rhoix}];
    fpdtype_t rl = ql[${rhoix}];

    // Compute the Roe-averaged enthalpy
    fpdtype_t H = (sqrtrl*(pr + ur[${Eix}])
                 + sqrtrr*(pl + ul[${Eix}]))
                / (sqrtrl*rr + sqrtrr*rl);

    // Roe average sound speed
    fpdtype_t u = (sqrtrl*vl[0] + sqrtrr*vr[0]) /
                  (sqrtrl + sqrtrr);

    fpdtype_t gamma = (sqrtrl*qhl[0] + sqrtrr*qhr[0]) /
                      (sqrtrl + sqrtrr);
    fpdtype_t a = sqrt((gamma - 1)*(H - 0.5*u*u));

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = u - a;
    fpdtype_t sr = u + a;
    fpdtype_t sstar = (pr - pl + rl*vl[0]*(sl - vl[0])
                               - rr*vr[0]*(sr - vr[0])) /
                      (rl*(sl - vl[0]) - rr*(sr - vr[0]));

    // Star state common factors
    fpdtype_t ul_com = (sl - vl[0]) / (sl - sstar);
    fpdtype_t ur_com = (sr - vr[0]) / (sr - sstar);

    // Star state mass
    fpdtype_t rusl = ul_com*rl;
    fpdtype_t rusr = ur_com*rr;
    // Star state species
% for n in range(ns):
    usl[${n}] = rusl*ql[${n}];
    usr[${n}] = rusr*qr[${n}];
% endfor

    // Star state momenetum
    usl[${vix}] = ul_com*rl*sstar;
    usr[${vix}] = ur_com*rr*sstar;

% for i in range(1, ndims):
    usl[${i + vix}] = ul_com*ul[${i + vix}];
    usr[${i + vix}] = ur_com*ur[${i + vix}];
% endfor

    // Star state energy
    usl[${Eix}] = ul_com*(ul[${Eix}] + (sstar - vl[0])*
                                (rl*sstar + pl/(sl - vl[0])));
    usr[${Eix}] = ur_com*(ur[${Eix}] + (sstar - vr[0])*
                                (rr*sstar + pr/(sr - vr[0])));

    // Output
% for i in range(nvars):
    fsl = fl[${i}] + sl*(usl[${i}] - ul[${i}]);
    fsr = fr[${i}] + sr*(usr[${i}] - ur[${i}]);
    nf[${i}] = (0 <= sl) ? fl[${i}] : (sl <= 0 && 0 <= sstar) ? fsl :
               (sstar <= 0 && 0 <= sr) ? fsr : fr[${i}];
% endfor
</%pyfr:macro>

<%include file='pyfr.solvers.mceuler.kernels.rsolvers.rsolve1d'/>
