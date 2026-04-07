<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.multicomp.makoutil' name='mc'/>

<% ns, vix, Eix, rhoix, pix, Tix = mc.thermix(c['ns'], ndims) %>

<%pyfr:macro name='compute_intestar' params='u, q, qh, intestar'>

  // For cpg, inte is always positive, just return qh3
  intestar = qh[3];

</%pyfr:macro>
