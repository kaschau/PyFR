from pyfr.multicomp.eos.base import BaseEOS
import itertools as it
import numpy as np
from numpy.polynomial import polynomial as nppoly
from scipy.optimize import least_squares

class MonotonicPolynomialFitter:
    def __init__(self, degree, interval):
        self.result_degree = degree
        self.q = degree - 1
        self.is_even = self.q % 2 == 0
        self.is_odd = not self.is_even
        self.interval = interval  # (a, b) or (a, None)

        # Determine degrees of p1 and p2 based on interval type and fit q
        q = self.q
        _, b = interval

        if b is None:
            # Semi-compact: [a, ∞)
            self.interval_type = "semi-compact"
            self.scaled_interval = (0, None)
        else:
            # Compact: [a, b]
            self.interval_type = "compact"
            self.scaled_interval = (0, 1)

        if q % 2 == 0:
            self.K = int(q / 2)
            self.p1_degree = self.K
            self.p2_degree = self.K - 1
        else:
            self.K = int((q - 1) / 2)
            self.p1_degree = self.p2_degree = self.K

        self.n_params = 1 + (self.p1_degree + 1) + (self.p2_degree + 1)

        # Scaling parameters (to be set during fit)
        self.x_scale = None
        self.x_offset = None
        self.y_scale = None
        self.y_offset = None

    def extract_params(self, theta):
        d1 = self.p1_degree + 1
        delta = theta[0]
        beta1 = np.array(theta[1 : 1 + d1])
        beta2 = np.array(theta[1 + d1 :])

        return delta, beta1, beta2

    def theta_to_coeffs(self, theta):
        """
        Convert parameter vector theta = [delta, beta1..., beta2...] to standard polynomial coefficients.

        Returns:
        beta: List of coefficients for the integrated monotonic polynomial.
              Represents p(x) = delta + 1.0 * ∫₀ˣ check_p(t) dt
        """
        delta, beta1, beta2 = self.extract_params(theta)

        # Build p1(x)^2 and p2(x)^2 via convolution
        p1_sq = nppoly.polymul(beta1, beta1)
        p2_sq = nppoly.polymul(beta2, beta2)

        # Determine weight function W(x) based on interval
        a, b = self.interval

        if self.interval_type == "semi-compact":
            # Multiply p2_sq by (x - a)
            p2_w = nppoly.polymul(p2_sq, [-a, 1])
            check_p = nppoly.polyadd(p1_sq, p2_w)
        else:
            if self.is_even:
                # W(x) = (x - a)(b - x) = -x^2 + (a + b)x - ab
                # So multiply p2_sq by a quadratic
                p2_w = nppoly.polymul(p2_sq, [-a * b, a + b, -1])
                check_p = nppoly.polyadd(p1_sq, p2_w)
            else:
                # p̌(x) = (x - a) p1(x)^2 + (b - x) p2(x)^2
                p1_w = nppoly.polymul(p1_sq, [-a, 1])
                p2_w = nppoly.polymul(p2_sq, [b, -1])
                check_p = nppoly.polyadd(p1_w, p2_w)

        # Integrate check_p(x) term-wise to get polynomial coefficients
        poly_coeffs = [delta]
        for i, c in enumerate(check_p):
            poly_coeffs.append(c / (i + 1))

        return poly_coeffs

    def evaluate_from_theta(self, x, theta):
        """
        Evaluate the monotonic polynomial at given x values.

        Parameters:
        x (float or array-like): Input value(s) at which to evaluate the polynomial.
        theta (list or array): Parameter vector [delta, beta1..., beta2...].

        Returns:
        y (float or array): Evaluated polynomial value(s).
        """
        coeffs = self.theta_to_coeffs(theta)
        return nppoly.polyval(x, coeffs)

    def residuals(self, theta, x, y):
        y_pred = self.evaluate_from_theta(x, theta)
        residuals = y - y_pred
        return residuals

    def gradient_rss(self, theta, x, y):
        """
        Compute gradient of RSS with respect to theta.

        Parameters:
        x (array-like): Input x values.
        y (array-like): Observed outputs.
        theta (list or array): [delta, beta1..., beta2...]

        Returns:
        grad (np.ndarray): Gradient vector of same shape as theta
        """

        _, beta1, beta2 = self.extract_params(theta)

        # residuals
        resid = self.residuals(theta, x, y)

        # Initialize gradient vector
        grad = np.zeros_like(theta)

        # ∂p/∂delta = 1
        grad[0] = -2 * np.sum(resid)

        # Build ∂p/∂βj1 terms
        grad_idx = 1
        for j in range(len(beta1)):
            dpdBj1 = np.zeros_like(x, dtype=np.float64)
            for k, Bk in enumerate(beta1):
                dpdBj1 += Bk * x ** (k + j + 1) / (k + j + 1)
            dpdBj1 *= 2

            grad[grad_idx] = -2 * np.sum(resid * dpdBj1)
            grad_idx += 1

        # Build ∂p/∂βj2 terms
        for j in range(len(beta2)):
            dpdBj2 = np.zeros_like(x, dtype=np.float64)
            for k, Bk in enumerate(beta2):
                dpdBj2 += Bk * x ** (k + j + 1) / (k + j + 1)
            dpdBj2 *= 2

            grad[grad_idx] = -2 * np.sum(resid * dpdBj2)
            grad_idx += 1

        return grad

    def initialize_theta(self, x, y):
        """
        Initialize theta0 based on simplified logic in the thesis (Section 2.5.2).
        No rescaling is done.
        Returns:
        - theta0: initial parameter vector [delta, beta1..., beta2...]
        """
        beta1 = np.zeros(self.p1_degree + 1)
        beta2 = np.zeros(self.p2_degree + 1)

        # Fit initial polynomial based on interval type and degree
        if self.interval_type == "semi-compact":
            if self.is_odd:
                # Fit y = A + Bx + Cx³
                X = np.column_stack([np.ones_like(x), x, x**3])
                coeffs = np.linalg.lstsq(X, y, rcond=None)[0]
                A, B, C = coeffs

                delta = A
                beta1[0] = np.sqrt(np.abs(B))
                if self.K >= 1:
                    beta1[1] = -np.sqrt(3 * np.abs(C))
                    beta2[0] = -np.sqrt(2 * np.abs(beta1[0] * beta1[1]))
            else:
                # Fit y = A + Bx + Cx²
                A, B, C = np.polynomial.Polynomial.fit(x, y, 2).convert().coef

                delta = A
                beta1[0] = np.sqrt(np.abs(B))
                beta2[0] = np.sqrt(2 * np.abs(C))

        elif self.interval_type == "compact":
            if self.is_odd:
                # Fit y = A + Bx + Cx³
                X = np.column_stack([np.ones_like(x), x, x**3])
                coeffs = np.linalg.lstsq(X, y, rcond=None)[0]
                A, B, C = coeffs

                delta = A
                beta1[0] = np.sqrt(np.abs(B))
                if self.K >= 1:
                    beta1[1] = -beta1[0] - np.sqrt(beta1[0] ** 2 + 3 * np.abs(C))
                    beta2[0] = np.sqrt(np.abs(beta1[1] ** 2 - 3 * np.abs(C)))
            else:
                # Fit y = A + Bx + Cx²
                A, B, C = np.polynomial.Polynomial.fit(x, y, 2).convert().coef

                delta = A
                beta2[0] = np.sqrt(np.abs(B))
                beta1[0] = -np.sqrt(np.abs(B) + 2 * np.abs(C))

        theta0 = np.concatenate([[delta], beta1, beta2])
        return theta0

    def scale_data(self, x, y):
        """
        Scale x and y data according to interval type.

        For semi-compact: x scaled so a=0 and max(x)=1
        For compact: x scaled so a=0 and b=1
        y scaled to range [-1, 1]

        Returns:
        x_scaled, y_scaled: scaled data arrays
        """
        a, b = self.interval

        # Scale x data
        if self.interval_type == "semi-compact":
            self.x_offset = a
            self.x_scale = np.max(x) - a
        else:  # compact
            self.x_offset = a
            self.x_scale = b - a
        x_scaled = (x - self.x_offset) / self.x_scale

        # Scale y data to [-1, 1]
        y_min, y_max = np.min(y), np.max(y)
        self.y_offset = (y_max + y_min) / 2
        self.y_scale = (y_max - y_min) / 2
        if self.y_scale == 0:
            self.y_scale = 1  # Avoid division by zero for constant data
        y_scaled = (y - self.y_offset) / self.y_scale

        return x_scaled, y_scaled

    def unscale_coefficients(self, coeffs_scaled):
        """
        Transform polynomial coefficients from scaled coordinates back to original coordinates.

        If p_scaled(x_scaled) = sum(c_i * x_scaled^i), then we need coefficients for
        p_original(x_original) such that p_original(x_original) = p_scaled(x_scaled) * y_scale + y_offset
        where x_scaled = (x_original - x_offset) / x_scale

        Parameters:
        coeffs_scaled: coefficients for polynomial in scaled coordinates

        Returns:
        coeffs_original: coefficients for polynomial in original coordinates
        """
        coeffs_scaled = np.array(coeffs_scaled)
        n = len(coeffs_scaled)
        coeffs_original = np.zeros(n)

        # Transform each coefficient
        for i in range(n):
            # Coefficient of x^i in original polynomial
            coeff_sum = 0
            for j in range(i, n):
                # Binomial expansion of ((x - x_offset) / x_scale)^j
                binomial_coeff = 1
                for k in range(i):
                    binomial_coeff *= (j - k) / (k + 1)

                term = (coeffs_scaled[j] * binomial_coeff *
                       ((-self.x_offset) ** (j - i)) / (self.x_scale ** j))
                coeff_sum += term

            coeffs_original[i] = coeff_sum * self.y_scale

        # Add y_offset to constant term
        coeffs_original[0] += self.y_offset

        return coeffs_original

    def max_relative_error(self, x, y, coeffs):
        """
        Calculate the maximum relative error over the fit data for unscaled coefficients.

        Parameters:
        x (array-like): Original x data used for fitting
        y (array-like): Original y data used for fitting
        coeffs (array-like): Unscaled polynomial coefficients

        Returns:
        float: Maximum relative error as |y_pred - y_actual| / |y_actual|
        """
        # Evaluate polynomial at original x points
        y_pred = nppoly.polyval(x, coeffs)

        # Calculate relative errors, avoiding division by zero
        relative_errors = np.abs(y_pred - y) / np.maximum(np.abs(y), 1e-12)

        return np.max(relative_errors)

    def is_monotonic(self, x, coeffs, is_scaled=True):
        """
        Verify that the fitted polynomial is indeed monotonic over the interval.
        """
        npts = 1000
        if is_scaled:
            interval = self.scaled_interval
        else:
            interval = self.interval

        rel = np.max(np.abs(x))
        if self.interval_type == 'semi-compact':
            x_min = interval[0]
            x_max = x.max() + rel
        elif self.interval_type == 'compact':
            x_min = interval[0]
            x_max = interval[1]
        x = np.linspace(x_min, x_max, npts)

        poly = nppoly.Polynomial(coeffs)
        polyder = poly.deriv()

        yder = polyder(x)

        if is_scaled:
            # We're in scaled coordinates
            scaled_pos_threshold = (1e-6 - self.y_offset) / self.y_scale
        else:
            # We're in original coordinates
            scaled_pos_threshold = 1e-6

        return np.all(yder > -1e-10), np.all(poly(x) > scaled_pos_threshold)

    def fit(self, x, y):
        # Scale the input data
        x_scaled, y_scaled = self.scale_data(x, y)

        # first test unconstrainted and see
        poly = nppoly.Polynomial.fit(x_scaled, y_scaled, self.result_degree)
        coef_scaled = poly.convert().coef
        mono, pos = self.is_monotonic(x_scaled, coef_scaled, is_scaled=True)
        if mono and pos:
            # Unscale coefficients and return
            return self.unscale_coefficients(coef_scaled)

        # otherwise we need to solve
        theta0 = self.initialize_theta(x_scaled, y_scaled)

        def rss(theta):
            return np.sum(self.residuals(theta, x_scaled, y_scaled)**2)

        def jac_func(theta):
            return self.gradient_rss(theta, x_scaled, y_scaled)

        lb = np.full(self.n_params, -np.inf)
        # Scale the lower bound for theta[0] to scaled coordinates
        # Original constraint: p(x >= 0) > 0 (small positive value)
        # In scaled coords: p_scaled >= (constraint - y_offset) / y_scale
        scaled_lower_bound = (1e-6 - self.y_offset) / self.y_scale
        lb[0] = scaled_lower_bound
        ub = np.full(self.n_params, np.inf)
        bounds = (lb, ub)
        result = least_squares(
            rss,
            theta0,
            jac=jac_func,
            gtol=1e-5,
            max_nfev=10000,
            bounds=bounds
        )

        theta_res = result.x
        coef_scaled = self.theta_to_coeffs(theta_res)

        # Unscale coefficients and return
        return self.unscale_coefficients(coef_scaled)


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
                # Fit cp
                m = np.where(Ts <= N7[n, 0], 8, 1)
                cp_ref = (       N7[n, m + 0]
                       + Ts*(N7[n, m + 1]
                       + Ts*(N7[n, m + 2]
                       + Ts*(N7[n, m + 3]
                       + Ts*(N7[n, m + 4] )))))

                # Adaptive monotonic polynomial fitting
                best_error = np.inf
                for degree in range(5):
                    fitter = MonotonicPolynomialFitter(degree, (0, None))
                    coeffs = list(fitter.fit(Ts, cp_ref))
                    error = fitter.max_relative_error(Ts, cp_ref, coeffs)
                    import matplotlib.pyplot as plt
                    plt.plot(Ts, cp_ref, label="ref")
                    cp_new = nppoly.polyval(Ts, coeffs)
                    plt.plot(Ts, cp_new, '--', label="new")
                    plt.title(f'c_p {consts['names'][n]}')
                    plt.legend()
                    plt.show()
                    if error < 0.01:
                        best_coeffs = coeffs
                        break
                    elif error < best_error:
                        best_coeffs = coeffs
                        best_error = error
                coeffs = best_coeffs

                h_ref  = (  Ts*(N7[n, m + 0]
                          + Ts*(N7[n, m + 1] / 2.0
                          + Ts*(N7[n, m + 2] / 3.0
                          + Ts*(N7[n, m + 3] / 4.0
                          + Ts*(N7[n, m + 4] / 5.0))))) + N7[n, m + 5])

                # Enthalpy from monotonic polynomial (integrated analytically)
                h_new = np.zeros_like(Ts)
                for i,coef in enumerate(coeffs):
                    h_new += coef * Ts**(i+1) / (i+1)

                # Find integration constant to match reference enthalpy
                a5 = np.mean(h_ref - h_new)
                coeffs.append(a5)

                # entropy integration constant
                s_ref = (  np.log(Ts)*N7[n, m + 0]
                                +(Ts *(N7[n, m + 1]
                                + Ts *(N7[n, m + 2] / 2.0
                                + Ts *(N7[n, m + 3] / 3.0
                                + Ts *(N7[n, m + 4] / 4.0))))) + N7[n, m + 6])

                # Entropy from monotonic polynomial
                s_new = coeffs[0] * np.log(Ts)
                for i in range(1, len(coeffs)-1):
                    s_new += coeffs[i] * Ts**i / i
                # Find integration constant to match reference entropy
                a6 = np.mean(s_ref - s_new)
                coeffs.append(a6)

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

        p, T = np.asarray(pris[0]), np.asarray(pris[ndims + 1])

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
            else:
                # Strict mode: use temperature-dependent NASA polynomials
                m = T <= consts['T_cutoff'][n]
                h_species = np.empty(T.shape)
                for idx, lh in zip([m, ~m], ['low','high']):
                    coeffs = consts[f'NASA7_T{lh}'][n]
                    # NASA polynomial enthalpy calculation
                    h_species[idx] = (  T[idx]*(coeffs[0]
                          + T[idx]*(coeffs[1] / 2.0
                          + T[idx]*(coeffs[2] / 3.0
                          + T[idx]*(coeffs[3] / 4.0
                          + T[idx]*(coeffs[4] / 5.0))))) + coeffs[5])

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

                    # Enthalpy polynomial: integrated C_p
                    h_species = 0.0
                    for i in range(len(coeffs) - 2):
                        h_species += coeffs[i] * T**(i+1) / (i+1)
                    h_species += coeffs[-2]  # Add enthalpy integration constant
                else:
                    # Strict mode: use temperature-dependent NASA polynomials
                    m = T <= consts['T_cutoff'][n]
                    cp_species = np.empty(T.shape)
                    h_species = np.empty(T.shape)
                    for idx, lh in zip([m, ~m], ['low','high']):
                        coeffs = consts[f'NASA7_T{lh}'][n]
                        # C_p calculation
                        cp_species[idx] = (     coeffs[0]
                                      + T[idx]*(coeffs[1]
                                      + T[idx]*(coeffs[2]
                                      + T[idx]*(coeffs[3]
                                      + T[idx]*(coeffs[4] )))))

                        # Enthalpy calculation
                        h_species[idx] = (  T[idx]*(coeffs[0]
                                            + T[idx]*(coeffs[1] /2.0
                                            + T[idx]*(coeffs[2] /3.0
                                            + T[idx]*(coeffs[3] /4.0
                                            + T[idx]*(coeffs[4] /5.0)))))
                                            +    coeffs[5])

                cp_species *= Ru / MW[n]
                cp += cp_species * Y
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