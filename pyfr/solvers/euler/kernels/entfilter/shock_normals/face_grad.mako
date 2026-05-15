<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Local divergence-theorem shock normal from a chosen scalar field
// (`fg_field`: density, pressure, or entropy) evaluated at every fpt of
// this cell using the self solution only -- no cross-cell information.
//
//   grad_ref ≈ (1/V_ref) sum_face n_{face,ref} * integral_face phi ds
//            ≈ (1/V_ref) sum_face n_{face,ref} * sum_{q on face} w_q * phi_q
//
// where phi_q is the chosen scalar at fpt q, computed by interpolating the
// upt solution via m0.  All face geometry (face_ref_normals, cell_ref_volume,
// facefpts, fpts_wts) is per-etype constant and baked at template time.
//
// For a smooth cell the per-fpt phi varies little along each face and
// opposite faces' contributions roughly cancel, giving a small grad_ref.
// For a shock-cut cell, phi varies sharply along the faces the shock
// crosses, the cancellation breaks, and grad_ref points perpendicular to
// the shock front.
//
// Output is in reference space; the caller converts to physical via smats^T
// (correct for isotropic Jacobian).
//
// Caveats per field choice:
//   density:  always defined; safest numerically.
//   pressure: well-defined but blows up at rho-undershoot fpts where the
//             kinetic energy term explodes.
//   entropy:  pow(rho, -gamma) undefined for rho<=0 -- avoid in cells with
//             pathological density undershoots.

<%pyfr:macro name='shock_normal_face_grad' params='u, m0, n_ref, n_mag'>
    // Compute phi at every fpt
    fpdtype_t phi[${nfpts}];
    {
        fpdtype_t uf[${nvars}];
        for (int fidx = 0; fidx < ${nfpts}; fidx++)
        {
            % for vidx in range(nvars):
            uf[${vidx}] = ${pyfr.dot('m0[fidx][{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
            % endfor

            % if fg_field == 'density':
            phi[fidx] = uf[0];
            % elif fg_field == 'pressure':
            {
                fpdtype_t rho = uf[0];
                fpdtype_t E = uf[${nvars - 1}];
                fpdtype_t ke = 0.5*(${' + '.join(
                    f'uf[{k}]*uf[{k}]' for k in range(1, ndims + 1))})/rho;
                phi[fidx] = ${c['gamma'] - 1}*(E - ke);
            }
            % elif fg_field == 'entropy':
            {
                fpdtype_t rho = uf[0];
                fpdtype_t E = uf[${nvars - 1}];
                fpdtype_t ke = 0.5*(${' + '.join(
                    f'uf[{k}]*uf[{k}]' for k in range(1, ndims + 1))})/rho;
                fpdtype_t p = ${c['gamma'] - 1}*(E - ke);
                phi[fidx] = (rho > 0 && p > 0)
                            ? p*pow(1.0/rho, ${c['gamma']}) : 0.0;
            }
            % else:
            <% raise ValueError(f"Unknown fg_field: {fg_field!r}") %>
            % endif
        }
    }

    // Divergence theorem: sum_face n * sum_fpt w_q * phi_q / V_ref
    fpdtype_t grad[${ndims}];
    % for d in range(ndims):
    grad[${d}] = ${' + '.join(
        f'phi[{fpt_idx}]*({face_ref_normals[f][d]*fpts_wts[fpt_idx]/cell_ref_volume})'
        for f, fpt_indices in enumerate(facefpts)
        for fpt_idx in fpt_indices
        if face_ref_normals[f][d] != 0
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
