<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>

<%include file='pyfr.solvers.mceuler.kernels.bcs.sub-in-mdot'/>
<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
