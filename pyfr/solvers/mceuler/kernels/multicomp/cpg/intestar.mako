<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='compute_intestar' params='u, q, qh, intestar'>

  // For cpg, inte is always positive, just return qh3
  intestar = qh[3]/q[${rhoix}];

</%pyfr:macro>
