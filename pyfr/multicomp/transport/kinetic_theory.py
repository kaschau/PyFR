from pyfr.multicomp.transport.base import BaseTransport
from pyfr.multicomp.MM_Tables import delta, tstar22, omega22_table, tstar, astar_table
import numpy as np
from scipy import interpolate as intrp


class KineticTheory(BaseTransport):
    name = "kinetic-theory"

    def __init__(self, cfg):
        super().__init__(cfg)

        # required properties
        self.input_props = [
            "MW",
            "well",
            "diam",
            "dipole",
            "polarize",
            "zrot",
            "geometry",
        ]

    def compute_consts(self, props, consts):
        self.consts = consts
        ns = consts["ns"]
        Ru = consts["Ru"]
        avogadro = consts["avogadro"]
        kb = consts["kb"]
        eps0 = consts["epsilon0"]

        eos = self.cfg.get("multi-component", "eos")
        if eos == "cpg":
            cp0 = props["cp0"]

            def cp_R(cp0, species_idx, T, MW, **kwargs):
                return cp0 * T / (Ru * MW)

        elif eos in ["tpg", "cubic"]:
            cp0 = [None for n in range(ns)]

            def cp_R(cp0, species_idx, T, MW, mode='strict', **kwargs):
                if mode == 'fast':
                    # Fast mode: use fitted coefficients
                    fast_coeff = kwargs['fast_coeff']
                    coeffs = fast_coeff[species_idx]
                    # Exclude integration constants (last two)
                    return sum([coeffs[i] * T ** i for i in range(len(coeffs) - 2)])
                elif mode == 'strict':
                    # Strict mode: use temperature-dependent NASA polynomials
                    T_cutoff = kwargs['T_cutoff']
                    NASA7_Thigh = kwargs['NASA7_Thigh']
                    NASA7_Tlow = kwargs['NASA7_Tlow']

                    # Select appropriate coefficient set based on temperature
                    if T <= T_cutoff[species_idx]:
                        coeffs = NASA7_Tlow[species_idx]
                    else:
                        coeffs = NASA7_Thigh[species_idx]

                    # coeffs = [c0, c1, c2, c3, c4, h_const, s_const]
                    return sum([coeffs[i] * T ** i for i in range(5)])
                else:
                    raise ValueError(f"Unknown mode: {mode}")

        maxdeg = 4
        prop_calc = self.cfg.get("multi-component", "property-calc", "strict")
        if prop_calc == "strict":
            fit = np.polynomial.Polynomial.fit
        else:
            from pyfr.multicomp.eos.base import poly_reduce as fit

        # Maximum and minumum temperatures to generate poly'l
        # NOTE: These ranges vary by input file in Cantera. It seems to set the
        # minTemp and maxTemp based on the min/max ranges of the NASA7 poly'l
        # data. In testing, this seems to explain the errors we sometimes get
        # in thermodynamic testing against Cantera.
        # (i.e. takes error from 1% to 0.001%).
        # For now we just use sensible values here.
        Tmin = self.cfg.getfloat("multi-component", "T-min", 300.0)
        Tmax = self.cfg.getfloat("multi-component", "T-max", 3500.0)
        # Generate range of temperatures
        npts = 400
        Ts = np.linspace(Tmin, Tmax, npts)

        # Collision integral interpolations
        intrp_o22 = intrp.RectBivariateSpline(tstar22, delta, omega22_table, kx=5, ky=5)
        intrp_Astar = intrp.RectBivariateSpline(tstar, delta, astar_table, kx=5, ky=5)

        # Get molecular mass
        MW = props["MW"]
        mass = np.array([M / avogadro for M in MW])

        # epsilon
        well = props["well"]
        # sigma
        diam = props["diam"]
        # mu
        dipole = props["dipole"]

        # see if molecule is polar
        polar = dipole > 0.0

        # alpha
        polarize = props["polarize"]

        # z_rot
        zrot = props["zrot"]

        # determine rotational DOF
        geom = props["geometry"]
        rotDOF = []
        for g in geom:
            if g == "atom":
                rotDOF.append(0.0)
            elif g == "linear":
                rotDOF.append(1.0)
            elif g == "nonlinear":
                rotDOF.append(1.5)

        ##########################################
        # Collision Parameters (reduced stuff)
        ##########################################
        r_mass = np.zeros((ns, ns))
        r_well = np.zeros((ns, ns))
        r_diam = np.zeros((ns, ns))
        r_dipole = np.zeros((ns, ns))

        r_deltastar = np.zeros((ns, ns))

        for n in range(ns):
            for n2 in range(n, ns):
                # reduced mass
                r_mass[n, n2] = mass[n] * mass[n2] / (mass[n] + mass[n2])
                # spheriacl collision diameter
                r_diam[n, n2] = 0.5 * (diam[n] + diam[n2])
                # effective well depth
                r_well[n, n2] = np.sqrt(well[n] * well[n2])
                # effective dipole moment
                r_dipole[n, n2] = np.sqrt(dipole[n] * dipole[n2])

                # reduced dipole delta*
                r_deltastar[n, n2] = (
                    0.5
                    * r_dipole[n, n2] ** 2
                    / (4 * np.pi * eps0 * r_well[n, n2] * r_diam[n, n2] ** 3)
                )

                # Correct for polarity
                if polar[n] == polar[n2]:
                    f_well = 1.0
                    f_diam = 1.0
                else:
                    kp, knp = (n, n2) if polar[n] else (n2, n)
                    d3np = diam[knp] ** 3
                    d3p = diam[kp] ** 3
                    alphastar = polarize[knp] / d3np
                    dipolestar = r_dipole[kp, kp] / np.sqrt(
                        4 * np.pi * eps0 * d3p * well[kp]
                    )
                    xi = 1.0 + 0.25 * alphastar * dipolestar**2 * np.sqrt(
                        well[kp] / well[knp]
                    )
                    f_well = xi**2
                    f_diam = xi ** (-1 / 6)

                r_well[n, n2] *= f_well
                r_diam[n, n2] *= f_diam

                # properties are symmetric
                r_mass[n2, n] = r_mass[n, n2]
                r_diam[n2, n] = r_diam[n, n2]
                r_well[n2, n] = r_well[n, n2]
                r_dipole[n2, n] = r_dipole[n, n2]
                r_deltastar[n2, n] = r_deltastar[n, n2]

        ##########################################
        # Viscosities
        ##########################################
        visc = np.zeros((npts, ns))
        for i, T in enumerate(Ts):
            Tstar = T * kb / well
            omga22 = intrp_o22(Tstar, r_deltastar.diagonal(), grid=False)
            visc[i, :] = (
                (5.0 / 16.0)
                * np.sqrt(np.pi * mass * kb * T)
                / (np.pi * diam**2 * omga22)
            )

        ##########################################
        # Thermal Conductivity
        ##########################################
        # NOTE The EOS has an effect on the transport properties via
        # the calculation of cp so if you use tpg you will use NASA7
        # to help compute thermal conductivities, if you use cpg you
        # will use constant cp to compute kappa.
        cond = np.zeros((npts, ns))
        for i, T in enumerate(Ts):
            for n in range(ns):
                Tstar = kb * 298.0 / well[n]
                fz_298 = (
                    1.0
                    + np.pi**1.5 / np.sqrt(Tstar) * (0.5 + 1.0 / Tstar)
                    + (0.25 * np.pi**2 + 2) / Tstar
                )

                Tstar = T * kb / well[n]

                omga22 = intrp_o22(Tstar, r_deltastar[n, n], grid=False)
                Astar = intrp_Astar(Tstar, r_deltastar[n, n], grid=False)
                omga11 = omga22 / Astar

                # self diffusion coeff
                diffcoeff = (
                    3.0
                    / 16.0
                    * np.sqrt(2.0 * np.pi / r_mass[n, n])
                    * (kb * T) ** 1.5
                    / (np.pi * diam[n] ** 2 * omga11)
                )

                f_int = MW[n] / (Ru * T) * diffcoeff / visc[i, n]
                cv_rot = rotDOF[n]
                A_factor = 2.5 - f_int
                fz_tstar = (
                    1.0
                    + np.pi**1.5 / np.sqrt(Tstar) * (0.5 + 1.0 / Tstar)
                    + (0.25 * np.pi**2 + 2) / Tstar
                )
                B_factor = zrot[n] * fz_298 / fz_tstar + 2.0 / np.pi * (
                    5 / 3 * rotDOF[n] + f_int
                )
                c1 = 2.0 / np.pi * A_factor / B_factor
                # Determine mode and get coefficient data
                if 'fast_coeff' in consts:
                    cv_int = cp_R(cp0[n], n, T, MW[n], mode='fast', fast_coeff=consts['fast_coeff']) - 2.5 - cv_rot
                else:
                    cv_int = cp_R(cp0[n], n, T, MW[n], mode='strict',
                                 T_cutoff=consts['T_cutoff'],
                                 NASA7_Thigh=consts['NASA7_Thigh'],
                                 NASA7_Tlow=consts['NASA7_Tlow']) - 2.5 - cv_rot
                f_rot = f_int * (1.0 + c1)
                f_trans = 2.5 * (1.0 - c1 * cv_rot / 1.5)
                cond[i, n] = (
                    visc[i, n]
                    / MW[n]
                    * Ru
                    * (f_trans * 1.5 + f_rot * cv_rot + f_int * cv_int)
                )

        ##########################################
        # Binary Diffusion
        ##########################################
        diff = np.zeros((npts, ns, ns))
        for n in range(ns):
            for j in range(n, ns):
                for i, T in enumerate(Ts):
                    Tstar = T * kb / r_well[j, n]

                    omga22 = intrp_o22(Tstar, r_deltastar[j, n], grid=False)
                    Astar = intrp_Astar(Tstar, r_deltastar[j, n], grid=False)
                    omga11 = omga22 / Astar

                    # To get pressure dependence, we evaluate the coeff at
                    # unit pressure then when we actually NEED the coeff,
                    # we use divide by the real pressure
                    diffcoeff = (
                        3.0
                        / 16.0
                        * np.sqrt(2.0 * np.pi / r_mass[n, j])
                        * (kb * T) ** 1.5
                        / (np.pi * r_diam[j, n] ** 2 * omga11)
                    )

                    diff[i, n, j] = diffcoeff
                    diff[i, j, n] = diff[i, n, j]

        # Create and set the polynoial coefficients
        logTs = np.log(Ts)
        sqrtTs = np.sqrt(Ts)

        # We fit the visc pol'y to the sqrtT as visc is proportional to sqrtT
        # we also reverse the numpy poly'l so lowest order is first
        visc = np.sqrt(visc / sqrtTs[:, None])
        w = 1.0 / (visc**2)
        consts["muPoly"] = [
            list(fit(logTs, visc[:, k], maxdeg, w=w[:, k]).convert().coef) for k in range(ns)
        ]

        # We fit the cond pol'y to the sqrtT as cond is proportional to sqrtT
        # we also reverse the numpy poly'l so lowest order is first
        cond = cond / np.sqrt(Ts[:, None])
        w = 1.0 / (cond**2)
        consts["kappaPoly"] = [
            list(fit(logTs, cond[:, k], maxdeg, w=w[:, k]).convert().coef) for k in range(ns)
        ]

        Dij = []
        diff = diff / Ts[:, None, None] ** 1.5
        w = 1.0 / (diff**2)
        for n in range(ns):
            for n2 in range(n, ns):
                poly = fit(logTs, 1.0/diff[:, n, n2], maxdeg, w=w[:, n, n2]).convert().coef
                Dij.append(list(poly))

        consts["DijPoly"] = Dij

        # MW should already be populated from the eos, but we redo it here anyway
        consts["MW"] = MW
