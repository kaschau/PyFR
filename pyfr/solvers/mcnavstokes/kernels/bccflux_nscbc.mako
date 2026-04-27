<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${mcf.trans}'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.flux'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.nscbc'/>
<%namespace file='pyfr.solvers.mcnavstokes.kernels.nscbc' import='nscbc_body'/>

<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u_upts='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              u_fpts='inout view fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              gradu_upts='in view fpdtype_t[${str(ndims*nupts)}][${str(nvars)}]'
              smats_upts='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              jacs_ffpts='in fpdtype_t[${str(nfacefpts)}]'>
${nscbc_body()}
</%pyfr:kernel>
