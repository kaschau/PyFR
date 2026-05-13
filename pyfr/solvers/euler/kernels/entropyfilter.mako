<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Entry point: dispatch to the named cascade, which provides the
// entropyfilter kernel body.  See pyfr/solvers/euler/kernels/entfilter/.

<%include file='pyfr.solvers.euler.kernels.entfilter.cascades.${cascade}'/>
