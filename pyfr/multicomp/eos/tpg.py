from pyfr.multicomp.eos.base import BaseEOS
import itertools as it
import numpy as np
from scipy.optimize import minimize, LinearConstraint
from scipy.linalg import lstsq

def fit_monotonic_polynomial(x, y, degree):
    """
    Fit a monotonic polynomial to data using constrained optimization.

    This function fits a polynomial of specified odd degree that is guaranteed
    to be monotonically increasing over the data range.

    Parameters:
    -----------
    x : array-like
        Input x values
    y : array-like
        Input y values
    degree : int
        Odd degree of the polynomial

    Returns:
    --------
    coefficients : ndarray
        Polynomial coefficients [c0, c1, c2, ..., c_degree]
        where polynomial is c0 + c1*x + c2*x^2 + ... + c_degree*x^degree

    Raises:
    -------
    ValueError
        If degree is not odd or if x and y have different lengths

    Example:
    --------
    >>> x = np.array([1, 2, 3, 4, 5])
    >>> y = np.array([1, 2.1, 3.3, 4.8, 6.2])
    >>> coefs = fit_monotonic_polynomial(x, y, degree=3)
    >>> print("Coefficients:", coefs)
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if degree % 2 == 0:
        raise ValueError("Degree must be odd for monotonic polynomials")

    if len(x) != len(y):
        raise ValueError("x and y must have the same length")

    # Check for constant or approximately constant data
    y_range = np.max(y) - np.min(y)
    y_mean = np.mean(y)

    # If y values are constant or nearly constant (within 5% relative tolerance), return constant polynomial
    relative_tolerance = 0.05  # 5%
    absolute_tolerance = 1e-12  # For near-zero values

    if y_range <= max(relative_tolerance * abs(y_mean), absolute_tolerance):
        # Return constant polynomial: f(x) = mean(y)
        coefficients = np.zeros(degree + 1)
        coefficients[0] = y_mean
        return coefficients

    # Sort data by x values if not already sorted
    if not np.all(x[:-1] <= x[1:]):
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        y = y[sort_idx]

    n_coef = degree + 1

    # Create Vandermonde matrix for polynomial fitting
    A = np.vander(x, n_coef, increasing=True)

    # Set up monotonicity constraints
    # For monotonicity, we need the derivative to be non-negative everywhere
    # The derivative is: c1 + 2*c2*x + 3*c3*x^2 + ...

    # Reduce constraint points for speed (still effective for monotonicity)
    n_constraint_points = min(20, max(10, len(x)))  # Much fewer points
    x_constraint = np.linspace(x.min(), x.max(), n_constraint_points)

    # Create constraint matrix for derivative >= 0
    constraint_matrix = np.zeros((n_constraint_points, n_coef))

    for i, xi in enumerate(x_constraint):
        for j in range(1, n_coef):  # Skip constant term (j=0)
            constraint_matrix[i, j] = j * (xi ** (j-1))

    # Constraint: derivative >= small positive value (for strict monotonicity)
    lower_bounds = np.full(n_constraint_points, 1e-6)
    upper_bounds = np.full(n_constraint_points, np.inf)

    # Objective function: minimize sum of squared residuals
    def objective(coef):
        pred = np.dot(A, coef)
        return np.sum((y - pred) ** 2)

    # Gradient of objective function
    def objective_grad(coef):
        pred = np.dot(A, coef)
        return 2 * np.dot(A.T, pred - y)

    # Hessian of objective function (constant for linear least squares)
    def objective_hess(coef):
        return 2 * np.dot(A.T, A)

    # Linear constraint for monotonicity
    linear_constraint = LinearConstraint(constraint_matrix, lower_bounds, upper_bounds)

    # Fast path: Try numpy solutions first
    # For degree 1: Use fast numpy polyfit with monotonicity
    if degree == 1:
        return fit_monotonic_linear_numpy(x, y)

    # Fast path: Check if unconstrained solution is already monotonic
    coef_unconstrained, _, _, _ = lstsq(A, y)
    if check_monotonicity(coef_unconstrained, (x.min(), x.max())):
        return coef_unconstrained

    # Fast path: Try constrained linear first (often good enough)
    if degree >= 3:
        linear_coeffs = fit_monotonic_linear_numpy(x, y)
        y_linear = linear_coeffs[0] + linear_coeffs[1] * x
        linear_error = np.max(np.abs(y - y_linear)) / np.max(np.abs(y))

        # If linear fit is already pretty good, use it instead of higher degree
        if linear_error < 0.02:  # 2% tolerance
            padded_coeffs = np.zeros(degree + 1)
            padded_coeffs[:2] = linear_coeffs
            return padded_coeffs

    # Try smart initial guess based on data characteristics
    coef_init = get_smart_initial_guess(x, y, degree, coef_unconstrained)

    # Quick check: if smart guess is already good enough and monotonic
    if check_monotonicity(coef_init, (x.min(), x.max())):
        y_fit = np.dot(A, coef_init)
        error = np.max(np.abs(y - y_fit)) / np.max(np.abs(y))
        if error < 0.01:  # 1% error tolerance for quick acceptance
            return coef_init

    # Use faster optimization with looser tolerances
    result = minimize(objective, coef_init, method='SLSQP',
                     jac=objective_grad,
                     constraints={'type': 'ineq',
                                'fun': lambda coef: np.dot(constraint_matrix, coef) - lower_bounds,
                                'jac': lambda coef: constraint_matrix},
                     options={'maxiter': 50, 'ftol': 1e-6, 'disp': False})  # Much looser tolerances

    if not result.success:
        print(f"Warning: Optimization may not have converged: {result.message}")

    return result.x

def get_smart_initial_guess(x, y, degree, coef_unconstrained):
    """Generate smart initial guess for constrained optimization"""

    # Start with unconstrained solution
    coef_init = coef_unconstrained.copy()

    # Strategy 1: If unconstrained has negative derivative terms, set them to small positive
    for i in range(1, len(coef_init)):
        if coef_init[i] < 0:
            coef_init[i] = 1e-6

    # Strategy 2: For higher degrees, use numpy to get data trends quickly
    if degree >= 3:
        # Get overall trend using fast numpy polyfit
        y_trend = np.polyfit(x, y, 1)[0]  # Overall slope

        if y_trend > 0:
            # Data is increasing, bias towards positive linear term
            coef_init[1] = max(coef_init[1], 0.1 * y_trend)

        # Damp higher-order terms that might cause oscillations
        coef_init[2:] *= 0.1  # Vectorized operation

    # Strategy 3: Ensure reasonable magnitude
    y_scale = np.std(y)
    x_scale = np.std(x)

    for i in range(1, len(coef_init)):
        # Scale coefficients to reasonable range
        expected_scale = y_scale / (x_scale ** i)
        if abs(coef_init[i]) > 10 * expected_scale:
            coef_init[i] = np.sign(coef_init[i]) * expected_scale

    return coef_init

def check_monotonicity(coefficients, x_range, n_points=100):
    """
    Check if polynomial is monotonic over given range.

    Parameters:
    -----------
    coefficients : array-like
        Polynomial coefficients
    x_range : tuple
        (x_min, x_max) range to check
    n_points : int
        Number of points to check

    Returns:
    --------
    is_monotonic : bool
        True if polynomial is monotonic (derivative >= 0) over the range
    """
    x_test = np.linspace(x_range[0], x_range[1], n_points)

    # Compute derivative coefficients
    deriv_coefs = [i * coefficients[i] for i in range(1, len(coefficients))]

    if len(deriv_coefs) == 0:
        return True  # Constant function is monotonic

    # Evaluate derivative
    deriv_values = np.zeros_like(x_test)
    for i, coef in enumerate(deriv_coefs):
        deriv_values += coef * (x_test ** i)

    # Check if derivative is non-negative (allowing small numerical errors)
    return np.all(deriv_values >= -1e-10)

class tpgEOS(BaseEOS):
    name = 'tpg'
    def __init__(self, cfg):
        super().__init__(cfg)

        self.input_props = [
            'MW',
            'NASA7',
        ]

    def compute_consts(self, props, consts):
        self.consts = consts
        consts['MW'] = props['MW']

        # Determine if we are using strict property cals or fast
        prop_calc = self.cfg.get('multi-component', 'property-calc', 'strict')
        if prop_calc == 'strict':
            # Restructure NASA7 data for cleaner access
            N7 = props['NASA7']
            ns = consts['ns']

            # Separate temperature cutoffs and coefficient arrays
            consts['T_cutoff'] = N7[:, 0]  # Temperature cutoff for each species

            # High temperature range coefficients (columns 1-7)
            consts['NASA7_Thigh'] = []
            for n in range(ns):
                coeffs = N7[n, 1:8].tolist()  # Convert to list for easier access
                consts['NASA7_Thigh'].append(coeffs)

            # Low temperature range coefficients (columns 8-14)
            consts['NASA7_Tlow'] = []
            for n in range(ns):
                coeffs = N7[n, 8:15].tolist()  # Convert to list for easier access
                consts['NASA7_Tlow'].append(coeffs)
        elif prop_calc == 'fast':
            N7 = props['NASA7']
            Tmin = self.cfg.getfloat('multi-component', 'T-min', 300.0)
            Tmax = self.cfg.getfloat('multi-component', 'T-max', 3500.0)
            Ts = np.linspace(Tmin, Tmax, 500)
            ns = consts['ns']

            # Store fitted coefficients as list of lists for cleaner access
            consts['fast_coeff'] = []

            for n in range(ns):
                # Determine which temperature range to use
                m = np.where(Ts <= N7[n, 0], 8, 1)

                w = np.ones(500)

                # Fit cp
                cp = (       N7[n, m + 0]
                       + Ts*(N7[n, m + 1]
                       + Ts*(N7[n, m + 2]
                       + Ts*(N7[n, m + 3]
                       + Ts*(N7[n, m + 4] )))))

                cp_poly = np.polynomial.Polynomial.fit(Ts, cp, 4, w=w)
                coeffs = list(cp_poly.convert().coef)

                # import matplotlib.pyplot as plt
                # plt.plot(Ts, cp, label="ref")
                # cp_new = (   coeffs[0]
                #        + Ts*(coeffs[1]
                #        + Ts*(coeffs[2]
                #        + Ts*(coeffs[3]
                #        + Ts*(coeffs[4] )))))
                # plt.plot(Ts, cp_new, '--', label="new")
                # plt.title(f'c_p {consts['names'][n]}')
                # plt.legend()
                # plt.show()

                h  = (  Ts*(N7[n, m + 0]
                      + Ts*(N7[n, m + 1] / 2.0
                      + Ts*(N7[n, m + 2] / 3.0
                      + Ts*(N7[n, m + 3] / 4.0
                      + Ts*(N7[n, m + 4] / 5.0))))) + N7[n, m + 5])
                h_new  = (  Ts*(coeffs[0]
                          + Ts*(coeffs[1] / 2.0
                          + Ts*(coeffs[2] / 3.0
                          + Ts*(coeffs[3] / 4.0
                          + Ts*(coeffs[4] / 5.0))))))

                # Use smallest h as integration constant
                hmin = np.argmin(np.abs(h))
                a5 = h[hmin] - h_new[hmin]
                coeffs.append(a5)

                # plt.plot(Ts, h, label="ref")
                # plt.plot(Ts, h_new + a5, "--", label="new")
                # plt.title(f'Enthalpy {consts['names'][n]}')
                # plt.legend()
                # plt.show()

                # entropy integration constant
                s = (  np.log(Ts)*N7[n, m + 0]
                            +(Ts *(N7[n, m + 1]
                            + Ts *(N7[n, m + 2] / 2.0
                            + Ts *(N7[n, m + 3] / 3.0
                            + Ts *(N7[n, m + 4] / 4.0))))) + N7[n, m + 6])

                s_new = (  np.log(Ts)*coeffs[0]
                               + (Ts*(coeffs[1]
                               +  Ts*(coeffs[2] / 2.0
                               +  Ts*(coeffs[3] / 3.0
                               +  Ts*(coeffs[4] / 4.0))))))
                smin = np.argmin(np.abs(s))
                a6 = s[smin] - s_new[smin]
                coeffs.append(a6)

                # plt.plot(Ts, s, label="ref")
                # plt.plot(Ts, s_new + a6, "--", label="new")
                # plt.title(f'Entropy {consts['names'][n]}')
                # plt.legend()
                # plt.show()

                # Store coefficients for this species
                consts['fast_coeff'].append(coeffs)
        else:
            raise ValueError(f'Unknown property-calc method "{prop_calc}".')


    def pri_to_con(self, pris):
        consts = self.consts
        Ru = consts['Ru']
        MW = consts['MW']
        ns = consts['ns']
        ndims = len(pris) - (ns - 1) - 2

        p, T = pris[0], pris[ndims + 1]

        # Compute ns species
        Yns = 1.0 - sum(pris[ndims+2::])

        # Check mass fractions all 0<Y<1
        self.validate_Y_ics(it.chain(pris[ndims+2::],[Yns]))

        # Compute mixture properties
        Rmix = 0.0
        for n, Y in enumerate(it.chain(pris[ndims+2::],[Yns])):
            Rmix += Y/consts['MW'][n]
        Rmix *= consts['Ru']

        # Compute h
        h = 0.0
        for n, Y in enumerate(it.chain(pris[ndims+2::],[Yns])):
            if 'fast_coeff' in consts:
                # Fast mode: use fitted coefficients
                coeffs = consts['fast_coeff'][n]
                # Integrate polynomial: c0*T + c1*T^2/2 + c2*T^3/3 + ... + h_const
                h_species = 0.0
                for i in range(len(coeffs) - 2):  # Exclude integration constants
                    h_species += coeffs[i] * T**(i+1) / (i+1)
                h_species += coeffs[-2]  # Add enthalpy integration constant
                h_species *= Ru / MW[n]
                h += h_species * Y
            else:
                # Strict mode: use temperature-dependent NASA polynomials
                if T <= consts['T_cutoff'][n]:
                    coeffs = consts['NASA7_Tlow'][n]
                else:
                    coeffs = consts['NASA7_Thigh'][n]

                # NASA polynomial enthalpy calculation
                h_species = (  T*(coeffs[0]
                      + T*(coeffs[1] / 2.0
                      + T*(coeffs[2] / 3.0
                      + T*(coeffs[3] / 4.0
                      + T*(coeffs[4] / 5.0))))) + coeffs[5])
                h_species *= Ru / MW[n]
                h += h_species * Y

        # Compute density
        rho = p/(Rmix*T)

        # Multiply velocity components by rho
        rhovs = [rho * c for c in pris[1 : ndims + 1]]

        # Compute the total energy
        rhok = 0.5 * rho * sum(c * c for c in pris[1 : ndims + 1])
        rhoe = rho * h - p
        rhoE = rhoe + rhok

        # Species mass
        rhoYk = [rho * c for c in it.chain(pris[ndims+2::],[Yns])]

        return [*rhoYk, *rhovs, rhoE]

    def con_to_pri(self, cons):
        consts = self.consts
        Ru = consts['Ru']
        MW = consts['MW']
        ns = consts['ns']
        ndims = len(cons)-(ns-1)-2

        rhoY = cons[0 : ns]
        rho = sum(rhoY)
        rhoE = cons[-1]
        # Divide momentum components by rho
        vs = [rhov / rho for rhov in cons[ns : ns + ndims]]

        # Species Mass Fraction
        Yk = [rhoYk / rho for rhoYk in rhoY]

        # Compute mixture properties
        Rmix = 0.0
        for n, Y in enumerate(Yk):
            Rmix += Y / MW[n]
        Rmix *= Ru

        # Internal energu
        e = rhoE/rho - 0.5 * sum(v * v for v in vs)

        # Iterate on T, start at 300K
        T = np.ones(e.shape)*300.0
        for _ in range(10):
            h = 0.0
            cp = 0.0
            for n, Y in enumerate(Yk):
                if 'fast_coeff' in consts:
                    # Fast mode: use fitted coefficients
                    coeffs = consts['fast_coeff'][n]

                    # C_p polynomial: c0 + c1*T + c2*T^2 + c3*T^3 + c4*T^4
                    cp_species = 0.0
                    for i in range(len(coeffs) - 2):  # Exclude integration constants
                        cp_species += coeffs[i] * T**i
                    cp_species *= Ru / MW[n]
                    cp += cp_species * Y

                    # Enthalpy polynomial: integrated C_p
                    h_species = 0.0
                    for i in range(len(coeffs) - 2):
                        h_species += coeffs[i] * T**(i+1) / (i+1)
                    h_species += coeffs[-2]  # Add enthalpy integration constant
                    h_species *= Ru / MW[n]
                    h += h_species * Y
                else:
                    # Strict mode: use temperature-dependent NASA polynomials
                    if T <= consts['T_cutoff'][n]:
                        coeffs = consts['NASA7_Tlow'][n]
                    else:
                        coeffs = consts['NASA7_Thigh'][n]

                    # C_p calculation
                    cp_species = (     coeffs[0]
                           + T*(coeffs[1]
                           + T*(coeffs[2]
                           + T*(coeffs[3]
                           + T*(coeffs[4] )))))
                    cp_species *= Ru / MW[n]
                    cp += cp_species * Y

                    # Enthalpy calculation
                    h_species = (  T*(coeffs[0]
                          + T*(coeffs[1] /2.0
                          + T*(coeffs[2] /3.0
                          + T*(coeffs[3] /4.0
                          + T*(coeffs[4] /5.0)))))
                          +    coeffs[5])
                    h_species *= Ru / MW[n]
                    h += h_species * Y
            error = e - (h - Rmix * T)
            # Newtons Method
            T -= error / (-cp + Rmix)

        p = rho*Rmix*T

        return [p, *vs, T, *Yk[0:-1]]

    def diff_con_to_pri(self, cons, diff_cons):
        consts = self.consts
        Ru = consts['Ru']
        MW = consts['MW']
        ns = consts['ns']
        ndims = len(cons) - (ns - 1) - 2

        rhoYk = cons[0:ns]
        *rhouvw, rhoE = cons[ns::]
        diff_rhoY= diff_cons[0:ns]
        *diff_rhouvw, diff_rhoE = diff_cons[ns::]

        rho = sum(rhoYk)
        diff_rho = sum(diff_rhoY)

        # Compute primiatives
        pris = self.con_to_pri(cons)
        p = pris[0]
        uvw = pris[1:ndims+1]
        T = pris[ndims+1]

        # Divide rhoY by ρ
        Yk = [rhoY / rho for rhoY in rhoYk]

        # Compute the temperature, pressure
        e = rhoE / rho - 0.5 * sum(v * v for v in uvw)

        # Velocity gradients: ∂u⃗ = 1/ρ·[∂(ρu⃗) - u⃗·∂ρ]
        diff_uvw = [(diff_rhov - v*diff_rho) / rho
                    for diff_rhov, v in zip(diff_rhouvw, uvw)]

        # Species gradients: ∂Y⃗ = 1/ρ·[∂(ρY⃗) - Y⃗·∂ρ]
        diff_Yk = [(diff_rhoY - Y*diff_rho) / rho
                    for diff_rhoY, Y in zip(diff_rhoY, Yk)]

        # Begin building temperature gradient
        diff_T = 1.0/rho*(diff_rhoE - rhoE/rho*diff_rho) - sum([i*j for i,j in zip(uvw,diff_uvw)])

        Rmix = 0.0
        cp = 0.0
        for n, (Y, diff_Y) in enumerate(zip(Yk, diff_Yk)):
            Rmix += Y / MW[n]

            if 'fast_coeff' in consts:
                # Fast mode: use fitted coefficients
                coeffs = consts['fast_coeff'][n]

                # C_p polynomial
                cp_species = 0.0
                for i in range(len(coeffs) - 2):
                    cp_species += coeffs[i] * T**i
                cp_species *= Ru / MW[n]
                cp += cp_species * Y

                # Enthalpy calculation for energy balance
                hk = 0.0
                for i in range(len(coeffs) - 2):
                    hk += coeffs[i] * T**(i+1) / (i+1)
                hk += coeffs[-2]  # Add enthalpy integration constant
                hk *= Ru / MW[n]

                e_Y = hk - T*Ru/MW[n]
                diff_T -= e_Y*diff_Y
            else:
                # Strict mode: use temperature-dependent NASA polynomials
                if T <= consts['T_cutoff'][n]:
                    coeffs = consts['NASA7_Tlow'][n]
                else:
                    coeffs = consts['NASA7_Thigh'][n]

                # C_p calculation
                cp_species = (     coeffs[0]
                       + T*(coeffs[1]
                       + T*(coeffs[2]
                       + T*(coeffs[3]
                       + T*(coeffs[4] )))))
                cp_species *= Ru / MW[n]
                cp += cp_species * Y

                # Enthalpy calculation
                hk = (  T*(coeffs[0]
                      + T*(coeffs[1] /2.0
                      + T*(coeffs[2] /3.0
                      + T*(coeffs[3] /4.0
                      + T*(coeffs[4] /5.0)))))
                      +    coeffs[5])
                hk *= Ru / MW[n]

                e_Y = hk - T*Ru/MW[n]
                diff_T -= e_Y*diff_Y

        Rmix *= Ru
        diff_T /= (cp - Rmix)

        # Build pressure gradient
        diff_p = Rmix*T*diff_rho + rho*Rmix*diff_T + rho*T*Ru*sum([dY/M for dY,M in zip(diff_Yk,MW)])

        return [diff_p, *diff_uvw, diff_T, *diff_Yk[0:-1]]