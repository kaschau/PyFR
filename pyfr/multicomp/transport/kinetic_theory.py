from pyfr.multicomp.transport.base_transport import BaseTransport
from pyfr.multicomp.MM_Tables import (delta,
                                      tstar22,
                                      omega22_table,
                                      tstar,
                                      astar_table)
import numpy as np
from scipy import interpolate as intrp



class KineticTheory(BaseTransport):
    name = 'kinetic-theory'

    def __init__(self):
        super().__init__()

        self.input_props = {
            'MW': None,
            'well': None,
            'diam': None,
            'dipole': None,
            'polarize': None,
            'zrot': None,
            'geometry': None,
        }

        self.consts = {
            'MW' : None,
            'muPoly' : None,
            'kappaPoly' : None,
            'DPoly' : None,
        }

    @staticmethod
    def compute_consts(props, consts, eos):
        ns = len(consts['names'])
        Ru = consts['Ru']
        avogadro = consts['avogadro']
        kb = consts['kb']
        eps0 = consts['epsilon0']

        if eos == "cpg":
            cp0 = props["cp0"]
            NASA7 = [None for n in range(ns)]

            def cp_R(cp0, poly, T, MW):
                return cp0 * T / (Ru * MW)

        elif eos in ["tpg", "cubic"]:
            NASA7 = props["NASA7"]
            cp0 = [None for n in range(ns)]

            def cp_R(cp0, poly, T, MW):
                if T <= poly[0]:
                    return sum([poly[i + 1 + 7] * T ** (i) for i in range(5)])
                else:
                    return sum([poly[i + 1] * T ** (i) for i in range(5)])

        deg = 4
        # Maximum and minumum temperatures to generate poly'l
        # NOTE: These ranges vary by input file in Cantera. It seems to set the
        # minTemp and maxTemp based on the min/max ranges of the NASA7 poly'l
        # data. In testing, this seems to explain the errors we sometimes get
        # in thermodynamic testing against Cantera.
        # (i.e. takes error from 1% to 0.001%).
        # For now we just use sensible values here.
        Tmin = 200
        Tmax = 5000
        # Generate range of temperatures
        npts = 50
        Ts = np.linspace(Tmin, Tmax, npts)

        # Collision integral interpolations
        intrp_o22 = intrp.RectBivariateSpline(tstar22, delta,
                                              omega22_table, kx=5, ky=5)
        intrp_Astar = intrp.RectBivariateSpline(tstar, delta,
                                                astar_table, kx=5, ky=5)

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

        for i in range(ns):
            for j in range(i, ns):
                # reduced mass
                r_mass[i, j] = mass[i] * mass[j] / (mass[i] + mass[j])
                # spheriacl collision diameter
                r_diam[i, j] = 0.5 * (diam[i] + diam[j])
                # effective well depth
                r_well[i, j] = np.sqrt(well[i] * well[j])
                # effective dipole moment
                r_dipole[i, j] = np.sqrt(dipole[i] * dipole[j])

                # reduced dipole delta*
                r_deltastar[i, j] = (
                    0.5
                    * r_dipole[i, j] ** 2
                    / (4 * np.pi * eps0 * r_well[i, j] * r_diam[i, j] ** 3)
                )

                # Correct for polarity
                if polar[i] == polar[j]:
                    f_well = 1.0
                    f_diam = 1.0
                else:
                    kp, knp = (i, j) if polar[i] else (j, i)
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

                r_well[i, j] *= f_well
                r_diam[i, j] *= f_diam

                # properties are symmetric
                r_mass[j, i] = r_mass[i, j]
                r_diam[j, i] = r_diam[i, j]
                r_well[j, i] = r_well[i, j]
                r_dipole[j, i] = r_dipole[i, j]
                r_deltastar[j, i] = r_deltastar[i, j]

        ##########################################
        # Viscosities
        ##########################################
        visc = np.zeros((npts, ns))
        for i, T in enumerate(Ts):
            Tstar = T * kb / well
            omga22 = intrp_o22(Tstar, r_deltastar.diagonal(), grid=False)
            visc[i, :] = (
                (5.0 / 16.0) * np.sqrt(np.pi * mass * kb * T)
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
            for k in range(ns):
                Tstar = kb * 298.0 / well[k]
                fz_298 = (
                    1.0
                    + np.pi**1.5 / np.sqrt(Tstar) * (0.5 + 1.0 / Tstar)
                    + (0.25 * np.pi**2 + 2) / Tstar
                )

                Tstar = T * kb / well[k]

                omga22 = intrp_o22(Tstar, r_deltastar[k, k], grid=False)
                Astar = intrp_Astar(Tstar, r_deltastar[k, k], grid=False)
                omga11 = omga22 / Astar

                # self diffusion coeff
                diffcoeff = (
                    3.0
                    / 16.0
                    * np.sqrt(2.0 * np.pi / r_mass[k, k])
                    * (kb * T) ** 1.5
                    / (np.pi * diam[k] ** 2 * omga11)
                )

                f_int = MW[k] / (Ru * T) * diffcoeff / visc[i, k]
                cv_rot = rotDOF[k]
                A_factor = 2.5 - f_int
                fz_tstar = (
                    1.0
                    + np.pi**1.5 / np.sqrt(Tstar) * (0.5 + 1.0 / Tstar)
                    + (0.25 * np.pi**2 + 2) / Tstar
                )
                B_factor = zrot[k] * fz_298 / fz_tstar + 2.0 / np.pi * (
                    5 / 3 * rotDOF[k] + f_int
                )
                c1 = 2.0 / np.pi * A_factor / B_factor
                cv_int = cp_R(cp0[k], NASA7[k], T, MW[k]) - 2.5 - cv_rot
                f_rot = f_int * (1.0 + c1)
                f_trans = 2.5 * (1.0 - c1 * cv_rot / 1.5)
                cond[i, k] = (
                    visc[i, k]
                    / MW[k]
                    * Ru
                    * (f_trans * 1.5 + f_rot * cv_rot + f_int * cv_int)
                )

        ##########################################
        # Binary Diffusion
        ##########################################
        diff = np.zeros((npts, ns, ns))
        for k in range(ns):
            for j in range(k, ns):
                for i, T in enumerate(Ts):
                    Tstar = T * kb / r_well[j, k]

                    omga22 = intrp_o22(Tstar, r_deltastar[j, k], grid=False)
                    Astar = intrp_Astar(Tstar, r_deltastar[j, k], grid=False)
                    omga11 = omga22 / Astar

                    # To get pressure dependence, we evaluate the coeff at
                    # unit pressure then when we actually NEED the coeff,
                    # we use divide by the real pressure
                    diffcoeff = (
                        3.0
                        / 16.0
                        * np.sqrt(2.0 * np.pi / r_mass[k, j])
                        * (kb * T) ** 1.5
                        / (np.pi * r_diam[j, k] ** 2 * omga11)
                    )

                    diff[i, k, j] = diffcoeff
                    diff[i, j, k] = diff[i, k, j]

        # Create and set the polynoial coefficients
        logTs = np.log(Ts)
        sqrtTs = np.sqrt(Ts)

        # We fit the visc pol'y to the sqrtT as visc is proportional to sqrtT
        # we also reverse the numpy poly'l so lowest order is first
        visc = np.sqrt(visc / sqrtTs[:, None])
        w = 1.0 / (visc**2)
        consts['muPoly'] = np.flip(
            np.array(
                [list(
                    np.polyfit(logTs, visc[:, k], deg=deg, w=w[:, k])
                    ) for k in range(ns)]
            ),
            -1,
        )

        # We fit the cond pol'y to the sqrtT as cond is proportional to sqrtT
        # we also reverse the numpy poly'l so lowest order is first
        cond = cond / np.sqrt(Ts[:, None])
        w = 1.0 / (cond**2)
        consts['kappaPoly'] = np.flip(
            np.array(
                [list(
                    np.polyfit(logTs, cond[:, k], deg=deg, w=w[:, k])
                    ) for k in range(ns)]
            ),
            -1,
        )

        Dij = np.empty((int((ns + 1)*ns/2), deg+1))
        diff = diff / Ts[:, None, None] ** 1.5
        w = 1.0 / (diff**2)
        icc = 0
        for n in range(ns):
            for n2 in range(n, ns):
                poly=np.polyfit(logTs, diff[:, n, n2], deg=deg, w=w[:, n, n2])
                Dij[icc,:] = np.flip(poly)
                icc += 1

        consts['DijPoly'] = Dij

        # MW should already be populated from the eos, but we redo it here anyway
        consts['MW'] = MW