<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='compute_intestar' params='u, q, qh, intestar'>

  // For cpg, inte is always positive, just return qh3
  intestar = qh[3];

</%pyfr:macro>
