<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>
<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='rsolve_1d' params='ul, ur, ql, qr, qhl, qhr, nf'>
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t usl[${nvars}], usr[${nvars}], fsl, fsr;
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t rhol = ul[0];
    fpdtype_t rhol = ur[0];
    fpdtype_t pl = ql[0];
    fpdtype_t pr = qr[0];

    ${pyfr.expand('inviscid_flux_1d', 'ul', 'fl', 'ql')};
    ${pyfr.expand('inviscid_flux_1d', 'ur', 'fr', 'qr')};

    // Wave speeds
    fpdtype_t cl = qhl[2];
    fpdtype_t cr = qhr[2];
    fpdtype_t sl = min(vl[0] - cl, vr[0] - cr);
    fpdtype_t sr = min(vl[0] + cl, vr[0] + cr);
    fpdtype_t sstar = (pr - pl + ul[${ndims+1}]*(sl - vl[0]) - ur[${ndims+1}]*(sr - vr[0])) /
                      (rhol*(sl - vl[0]) - rhor*(sr - vr[0]));

    fpdtype_t ul_com = (sl - vl[0]) / (sl - sstar);
    fpdtype_t ur_com = (sr - vr[0]) / (sr - sstar);

% for n in range(ns):
    usl[${Yix + n}] = ul_com*ul[${Yix + n}];
    usr[${Yix + n}] = ur_com*ur[${Yix + n}];
% endfor

    usl[0] = ul_com*rhol*sstar;
    usr[0] = ur_com*rhor*sstar;
% for i in range(1, ndims):
    usl[${i}] = ul_com*ul[${i}];
    usr[${i}] = ur_com*ur[${i}];
% endfor

    usl[${ndims+1}] = ul_com*(ul[${ndims+1}] + rhol*(sstar - vl[0])*(sstar + pl/(rhol*(sl - vl[0]))));
    usr[${ndims+1}] = ur_com*(ur[${ndims+1}] + rhor*(sstar - vr[0])*(sstar + pr/(rhor*(sr - vr[0]))));

    // Output
% for i in range(nvars):
    fsl = fl[${i}] + sl*(usl[${i}] - ul[${i}]);
    fsr = fr[${i}] + sr*(usr[${i}] - ur[${i}]);
    nf[${i}] = (0 <= sl) ? fl[${i}] : (sl <= 0 && 0 <= sstar) ? fsl : (sstar <= 0 && 0 <= sr) ? fsr : fr[${i}];
% endfor
</%pyfr:macro>

<%include file='pyfr.solvers.mpeuler.kernels.rsolvers.rsolve1d'/>