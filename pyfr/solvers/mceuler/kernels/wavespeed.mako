<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.multicomp.makoutil' name='mc'/>
<%include file='pyfr.solvers.baseadvec.kernels.smats'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<% smats = 'smats_l' if 'linear' in ktype else 'smats' %>
<% rcpdjac_v = 'rcpdjac_l' if 'linear' in ktype else 'rcpdjac' %>
<% ns, vix, Eix, rhoix, pix, Tix = mc.thermix(c['ns'], ndims) %>

<%pyfr:kernel name='wavespeed' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t'
              verts='in broadcast-col fpdtype_t[${str(nverts)}][${str(ndims)}]'
              upts='in broadcast-row fpdtype_t[${str(ndims)}]'
              wspd='out broadcast-col reduce(max) fpdtype_t'>
% if 'linear' in ktype:
    fpdtype_t ${smats}[${ndims}][${ndims}], djac;
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', smats, 'djac')};
    fpdtype_t ${rcpdjac_v} = 1/djac;
% endif

    // Compute thermodynamic state
    fpdtype_t q[${nvars + 2}];
    fpdtype_t qh[${4 + ns}];
    ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};

    fpdtype_t csnd = qh[2];

    fpdtype_t lam = 0;
% for i in range(ndims):
    lam += fabs(${' + '.join(f'({smats}[{i}][{j}]*{rcpdjac_v})*q[{vix + j}]'
                             for j in range(ndims))})
         + csnd*sqrt(${' + '.join(f'({smats}[{i}][{j}]*{rcpdjac_v})'
                                  f'*({smats}[{i}][{j}]*{rcpdjac_v})'
                                  for j in range(ndims))});
% endfor

    wspd = lam;
</%pyfr:kernel>
