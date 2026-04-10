<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-prims'/>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

%  for n,spn in enumerate(mcf.sp_names):
    qr[${n}] = ${c[spn]};
%  endfor

% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + vix}] = ${c[v]};
% endfor

    qr[${pix}] = ${c['p']};
    qr[${Tix}] = ${c['T']};

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>
