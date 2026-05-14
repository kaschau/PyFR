<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Reference-space divergence-theorem shock normal from per-face entmin.
//
//   grad_ref ≈ (1/V_ref) sum_f entmin_int[f] * n_{f,ref} * L_{f,ref}
//
// All face geometry is per-etype constants (face_ref_normals, face_ref_lengths,
// cell_ref_volume in tplargs), so the macro fully unrolls at template time
// with no per-element data needed beyond entmin_int.
//
// For axis-aligned linear elements this output is the (covariant) reference-
// space gradient direction, which equals the (contravariant) filter-ready
// `n_ref = J^-1 n_phys` direction up to renormalisation when J is isotropic
// (as it is for axis-aligned linear quads/hexes).  For general curvilinear
// meshes the limiter that eventually consumes this would need an additional
// J^-T J^-1 transform.  The TODO lives at the consumer, not here.

<%pyfr:macro name='shock_normal_face_entmin' params='entmin_int, n_ref, n_mag'>
    fpdtype_t grad[${ndims}];
    % for d in range(ndims):
    grad[${d}] = ${' + '.join(
        f'entmin_int[{f}]*({fn[d]*L/cell_ref_volume})'
        for f, (fn, L) in enumerate(zip(face_ref_normals, face_ref_lengths))
        if fn[d] != 0
    ) or '0.0'};
    % endfor

    fpdtype_t mag2 = ${' + '.join(f'grad[{d}]*grad[{d}]' for d in range(ndims))};
    n_mag = sqrt(mag2);

    if (n_mag > ${shock_normal_eps})
    {
        fpdtype_t inv_mag = 1.0/n_mag;
        % for d in range(ndims):
        n_ref[${d}] = grad[${d}]*inv_mag;
        % endfor
    }
    else
    {
        % for d in range(ndims):
        n_ref[${d}] = 0.0;
        % endfor
    }
</%pyfr:macro>
