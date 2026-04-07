<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.multicomp.makoutil' name='mc'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>

<% ns, vix, Eix, rhoix, pix, Tix = mc.thermix(c['ns'], ndims) %>

<%pyfr:macro name='rsolve' params='ul, ur, ql, qr, qhl, qhr, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr, nf_fl, nf_fr, nf_fsl, nf_fsr;
    fpdtype_t va[${ndims}];
    fpdtype_t usl[${nvars}], usr[${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'ql')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'qr')};

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

    // Get the normal left and right velocities
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};
    fpdtype_t al = qhl[2];
    fpdtype_t ar = qhr[2];
    // Compute the Roe-averaged velocity
    fpdtype_t nv = (sqrtrl*nvl + sqrtrr*nvr)
                 / (sqrtrl + sqrtrr);

    // Compute the Roe-averaged enthalpy
    fpdtype_t H = (sqrtrl*(pr + ur[${Eix}])
                 + sqrtrr*(pl + ul[${Eix}]))
                / (sqrtrl*rr + sqrtrr*rl);

    fpdtype_t inv_rar = 1 / (sqrtrl + sqrtrl);
% for i in range(ndims):
    va[${i}] = (vl[${i}]*sqrtrl + vr[${i}]*sqrtrr) * inv_rar;
% endfor
    fpdtype_t qq = ${pyfr.dot('va[{i}]', i=ndims)};

    // Roe average speed of sound
    fpdtype_t gamma = (sqrtrl*qhl[0] + sqrtrr*qhr[0]) /
                      (sqrtrl + sqrtrr);
    fpdtype_t a = sqrt((gamma - 1)*(H - 0.5*qq));

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = min(nv - a, nvl - al);
    fpdtype_t sr = max(nv + a, nvr + ar);
    fpdtype_t sstar = (pr - pl + rl*nvl*(sl - nvl)
                               - rr*nvr*(sr - nvr)) /
                      (rl*(sl - nvl) - rr*(sr - nvr));


    // Star state common factors
    fpdtype_t ul_com = (sl - nvl) / (sl - sstar);
    fpdtype_t ur_com = (sr - nvr) / (sr - sstar);

    // Star state mass
    fpdtype_t rusl = ul_com*rl;
    fpdtype_t rusr = ur_com*rr;
    // Star state species
% for n in range(ns):
    usl[${n}] = rusl*ql[${n}];
    usr[${n}] = rusr*qr[${n}];
% endfor

    // Star state momenetum
% for i in range(ndims):
    usl[${vix + i}] = rusl*(vl[${i}] + (sstar - nvl)*n[${i}]);
    usr[${vix + i}] = rusr*(vr[${i}] + (sstar - nvr)*n[${i}]);
%endfor

    // Star state energy
    usl[${Eix}] = ul_com*(ul[${Eix}] + (sstar - nvl)*
                                (rl*sstar + pl/(sl - nvl)));
    usr[${Eix}] = ur_com*(ur[${Eix}] + (sstar - nvr)*
                                (rr*sstar + pr/(sr - nvr)));

    // Output
% for i in range(nvars):
    nf_fl = ${' + '.join(f'n[{j}]*fl[{j}][{i}]' for j in range(ndims))};
    nf_fr = ${' + '.join(f'n[{j}]*fr[{j}][{i}]' for j in range(ndims))};
    nf_fsl = nf_fl + sl*(usl[${i}] - ul[${i}]);
    nf_fsr = nf_fr + sr*(usr[${i}] - ur[${i}]);
    nf[${i}] = (0 <= sl) ? nf_fl : (sl <= 0 && 0 <= sstar) ? nf_fsl :
               (sstar <= 0 && 0 <= sr) ? nf_fsr : nf_fr;
% endfor
</%pyfr:macro>