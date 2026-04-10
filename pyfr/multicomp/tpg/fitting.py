import numpy as np

from pyfr.multicomp.tpg.collision_tables import (
    delta, tstar22, omega22_table, astar_table
)


# -----------------------------------------------------------------------
# Polynomial fitting
# -----------------------------------------------------------------------

def fit_poly(x, y, tol=0, w=None, maxdeg=4):
    if tol > 0:
        for deg in range(maxdeg + 1):
            p = np.polynomial.Polynomial.fit(x, y, deg, w=w)
            coeffs = list(p.convert().coef)
            y_fit = np.polyval(coeffs[::-1], x)
            err = np.max(np.abs((y_fit - y) / y))
            if err < tol:
                return coeffs
        return coeffs
    else:
        p = np.polynomial.Polynomial.fit(x, y, maxdeg, w=w)
        return list(p.convert().coef)


# -----------------------------------------------------------------------
# NASA polynomial evaluation (Horner form)
# -----------------------------------------------------------------------

def eval_nasa7_cp(T, a):
    return a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*a[4])))


def eval_nasa7_h(T, a):
    return T*(a[0] + T*(a[1]/2 + T*(a[2]/3 + T*(a[3]/4 + T*a[4]/5)))) + a[5]


def eval_nasa7_s(T, a):
    return (a[0]*np.log(T)
            + T*(a[1] + T*(a[2]/2 + T*(a[3]/3 + T*a[4]/4))) + a[6])


def eval_nasa9_cp(T, a):
    T2 = T*T
    return a[0]/T2 + a[1]/T + a[2] + T*(a[3] + T*(a[4] + T*(a[5] + T*a[6])))


def eval_nasa9_h(T, a):
    return (-a[0]/T + a[1]*np.log(T) + T*(a[2] + T*(a[3]/2 + T*(a[4]/3
            + T*(a[5]/4 + T*a[6]/5)))) + a[7])


def eval_nasa9_s(T, a):
    T2 = T*T
    return (-a[0]/(2*T2) - a[1]/T + a[2]*np.log(T)
            + T*(a[3] + T*(a[4]/2 + T*(a[5]/3 + T*a[6]/4))) + a[8])


def eval_over_ranges(T, ranges, coeffs, eval_fn):
    result = np.empty_like(T)
    for j in range(len(coeffs)):
        if j == 0:
            mask = T <= ranges[j + 1]
        elif j == len(coeffs) - 1:
            mask = T > ranges[j]
        else:
            mask = (T > ranges[j]) & (T <= ranges[j + 1])
        if np.any(mask):
            result[mask] = eval_fn(T[mask], coeffs[j])
    return result


# -----------------------------------------------------------------------
# Collision integral interpolation (matches Cantera's quadInterp)
# -----------------------------------------------------------------------

_log_tstar22 = np.log(tstar22)


def _quad_interp(x0, x, y):
    dx21 = x[1] - x[0]
    dx32 = x[2] - x[1]
    dx31 = dx21 + dx32
    dy32 = y[2] - y[1]
    dy21 = y[1] - y[0]
    a = (dx21*dy32 - dy21*dx32) / (dx21*dx31*dx32)
    return a*(x0 - x[0])*(x0 - x[1]) + (dy21/dx21)*(x0 - x[1]) + y[1]


def _interp_table(ts_val, deltastar, xs_table, table, log_xs, offset=0):
    n = len(xs_table)

    i = 0
    for j in range(n):
        if ts_val < xs_table[j]:
            i = j
            break
    else:
        i = n

    i1 = max(i - 1, 0)
    i2 = i1 + 3
    if i2 > n:
        i2 = n
        i1 = i2 - 3

    vals = np.empty(3)
    for k, idx in enumerate(range(i1, i2)):
        if deltastar == 0.0:
            vals[k] = table[idx + offset, 0]
        else:
            ds = np.clip(deltastar, delta[0], delta[-1])
            p = np.polynomial.Polynomial.fit(delta, table[idx + offset, :],
                                              min(6, len(delta) - 1))
            vals[k] = p(ds)

    return _quad_interp(np.log(ts_val), log_xs[i1:i2], vals)


def omega22_at(ts_val, deltastar):
    return _interp_table(ts_val, deltastar, tstar22, omega22_table,
                         _log_tstar22)


def astar_at(ts_val, deltastar):
    return _interp_table(ts_val, deltastar, tstar22, astar_table,
                         _log_tstar22, offset=1)


def fit_collision_integrals(deltastar, tstar_min=None, tstar_max=None,
                            degree=8):
    """Pre-fit omega22 and A* to degree-8 polynomials in log(Tstar).

    Matches Cantera's fitCollisionIntegrals approach: interpolate the
    table in delta* at each Tstar row, trim to relevant range, then fit.

    Returns (o22_coeffs, astar_coeffs) as lists of polynomial coefficients.
    """
    def _fit_table(xs, table, offset=0):
        # Trim Tstar range
        if tstar_min is not None and tstar_max is not None:
            nmin = nmax = 0
            for i in range(len(xs)):
                if tstar_min > xs[i]:
                    nmin = i
                if tstar_max > xs[i]:
                    nmax = i
            nmax = min(nmax + 1, len(xs) - 1)
        else:
            nmin, nmax = 0, len(xs) - 1

        xs_trim = xs[nmin:nmax + 1]
        log_xs = np.log(xs_trim)

        vals = np.empty(len(xs_trim))
        for i in range(len(xs_trim)):
            idx = nmin + i
            if deltastar == 0.0:
                vals[i] = table[idx + offset, 0]
            else:
                ds = np.clip(deltastar, delta[0], delta[-1])
                p = np.polynomial.Polynomial.fit(
                    delta, table[idx + offset, :],
                    min(6, len(delta) - 1)
                )
                vals[i] = p(ds)

        deg = min(degree, len(xs_trim) - 1)
        p = np.polynomial.Polynomial.fit(log_xs, vals, deg)
        return list(p.convert().coef)

    o22_coeffs = _fit_table(tstar22, omega22_table)
    astar_coeffs = _fit_table(tstar22, astar_table, offset=1)
    return o22_coeffs, astar_coeffs


def eval_collision_poly(log_tstar, coeffs):
    """Evaluate a pre-fitted collision integral polynomial."""
    result = np.zeros_like(log_tstar)
    for i, c in enumerate(coeffs):
        result += c * log_tstar**i
    return result
