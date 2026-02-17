def thermix(context, ns, ndims):
    vix = ns
    Eix = ns + ndims
    rhoix = ns + ndims
    pix = rhoix + 1
    Tix = pix + 1
    return ns, vix, Eix, rhoix, pix, Tix

def nasa_cps(context, coeffs, Ru, MW):
    """
    Generate C_p polynomial expression from coefficient list.
    coeffs: list of polynomial coefficients [c0, c1, c2, c3, c4, ...]
    Last two coefficients are assumed to be h_const and s_const
    """
    # Use all coefficients except the last two (which are integration constants)
    poly_coeffs = coeffs[:-2]
    scaled_coeffs = [c * Ru / MW for c in poly_coeffs]

    # Build nested polynomial expression using Horner's method
    return f'({'+ T*('.join(str(c) for c in scaled_coeffs)+')'*(len(scaled_coeffs)-1)})'

def nasa_cpp(context, coeffs, Ru, MW):
    """
    Generate dC_p/dT polynomial expression from coefficient list.
    """
    # Derivative coefficients: multiply by power and shift down
    poly_coeffs = coeffs[:-2]  # Exclude integration constants
    deriv_coeffs = []
    for i in range(1, len(poly_coeffs)):
        deriv_coeffs.append(i * poly_coeffs[i] * Ru / MW)

    if len(deriv_coeffs) == 0:
        return '0.0'

    # Build nested polynomial expression
    return f'({'+ T*('.join(str(c) for c in deriv_coeffs)+')'*(len(deriv_coeffs)-1)})'

def nasa_hs(context, coeffs, Ru, MW):
    """
    Generate enthalpy polynomial expression from coefficient list.
    """
    poly_coeffs = coeffs[:-2]  # Exclude integration constants
    h_const = coeffs[-2]  # Second to last is enthalpy integration constant

    # Integrate C_p: divide by (i+1) for each term
    integrated_coeffs = []
    for i, c in enumerate(poly_coeffs):
        integrated_coeffs.append(c * Ru / MW / (i + 1))

    # Build nested polynomial expression
    poly_expr = f'T*({'+ T*('.join(str(c) for c in integrated_coeffs)+')'*len(integrated_coeffs)}'
    return f'({poly_expr} + {h_const * Ru/MW})'

def nasa_hi(context, coeffs):
    poly_coeffs = coeffs[:-2]  # Exclude integration constants
    h_const = coeffs[-2]  # Second to last is enthalpy integration constant
    integrated_coeffs = []
    for i, c in enumerate(poly_coeffs):
        integrated_coeffs.append(c / (i + 1))

    # Build nested polynomial expression
    poly_expr = f'({'+ T*('.join(str(c) for c in integrated_coeffs)+')'*len(integrated_coeffs)}'
    return f'({poly_expr} + {h_const}*Tinv)'

def nasa_gbs(context, coeffs):
    """
    Generate Gibbs free energy polynomial expression from coefficient list.
    G = H - T*S
    """
    poly_coeffs = coeffs[:-2]  # Exclude integration constants
    h_const = coeffs[-2]  # Enthalpy integration constant
    s_const = coeffs[-1]  # Entropy integration constant

    # G = H - T*S = (H/T) - S
    # Start with log term from entropy: c0*(1 - ln(T))
    log_coeff = poly_coeffs[0]
    terms = [f'{log_coeff} * (1.0 - logT)']

    # Add enthalpy constant term: h_const/T
    terms.append(f'{h_const} * Tinv')

    # Add polynomial terms: H/T - S integration
    if len(poly_coeffs) > 1:
        # For polynomial terms c_i*T^i:
        # H contribution: c_i*T^i/i -> c_i*T^(i-1)/i (divided by T)
        # S contribution: -c_i*T^i/i -> subtract this
        div_coeffs = []
        for i in range(1, len(poly_coeffs)):
            # Combined: c_i*T^(i-1)/i - c_i*T^i/i = c_i/i * (T^(i-1) - T^i) = -c_i*T^(i-1)*(T-1)/i
            # Simplified for Gibbs: just the integration pattern
            div_coeffs.append(poly_coeffs[i] / (i * (i + 1)))

        if div_coeffs:
            poly_expr = f'T*({'+ T*('.join(str(-c) for c in div_coeffs)+')'*len(div_coeffs)}'
            terms.append(poly_expr)

    # Subtract entropy constant
    terms.append(f'- {s_const}')

    return f'({" + ".join(terms).replace("+ -", "- ")})'

def nasa_s(context, coeffs, Ru, MW):
    """
    Generate entropy polynomial expression from coefficient list (with units).
    """
    poly_coeffs = coeffs[:-2]  # Exclude integration constants
    s_const = coeffs[-1]  # Last is entropy integration constant

    # Entropy: S = c0*ln(T) + c1*T + c2*T^2/2 + c3*T^3/3 + ...
    log_coeff = poly_coeffs[0] * Ru / MW

    if len(poly_coeffs) > 1:
        # Integration of polynomial terms (divide by power)
        integrated_coeffs = []
        for i in range(1, len(poly_coeffs)):
            integrated_coeffs.append(poly_coeffs[i] * Ru / MW / i)

        poly_expr = f'T*({'+ T*('.join(str(c) for c in integrated_coeffs)+')'*len(integrated_coeffs)}'
        return f'({log_coeff} * logT + {poly_expr} + {s_const * Ru/MW})'
    else:
        return f'({log_coeff} * logT + {s_const * Ru/MW})'
