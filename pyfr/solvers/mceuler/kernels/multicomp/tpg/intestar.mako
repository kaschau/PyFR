<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='compute_intestar' params='u, q, qh, intestar'>
    intestar = qh[3];
% for n in range(mcf.ns):
    intestar -= u[${n}] * ${mcf[n].h_ref};
% endfor
</%pyfr:macro>
