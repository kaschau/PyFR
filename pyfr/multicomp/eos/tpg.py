from pyfr.multicomp.eos.base import BaseEOS
import itertools as it
import numpy as np

def evaluate_polynomial(x, coefficients):
    x = np.asarray(x, dtype=float)
    result = np.zeros_like(x, dtype=float)

    for i, coef in enumerate(coefficients):
        result += coef * (x ** i)

    return result

def fit_adaptive_monotonic_polynomial(x, y, tolerance=0.01):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    for degree in [0,1,2,3,4]:
        # Fit unconstrained
        poly = np.polynomial.Polynomial.fit(x, y, degree)
        coeffs = list(poly.convert().coef)

        # Evaluate polynomial
        y_fit = evaluate_polynomial(x, coeffs)

        # Calculate L∞ norm (maximum absolute error)
        abs_error = np.max(np.abs(y_fit - y))
        rel_error = abs_error / np.max(np.abs(y)) if np.max(np.abs(y)) > 1e-12 else abs_error

        if rel_error < tolerance:
            return coeffs

    return coeffs


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
                coeffs = fit_adaptive_monotonic_polynomial(Ts, cp_ref)

                # import matplotlib.pyplot as plt
                # plt.plot(Ts, cp_ref, label="ref")
                # cp_new = evaluate_polynomial(Ts, coeffs)
                # plt.plot(Ts, cp_new, '--', label="new")
                # plt.title(f'c_p {consts['names'][n]}')
                # plt.legend()
                # plt.show()

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