from collections import namedtuple
from collections.abc import Iterable
from inspect import signature
import itertools as it
import re

from mako.runtime import capture, supports_caller
import numpy as np

import pyfr.nputil as nputil
import pyfr.util as util


def ndrange(context, *args):
    return util.ndrange(*args)


def ilog2range(context, x):
    return [2**i for i in range(x.bit_length() - 2, -1, -1)]


def npdtype_to_ctype(context, dtype):
    return nputil.npdtype_to_ctype(dtype)


def dot(context, a_, b_=None, /, **kwargs):
    ix, nd = util.first(kwargs.items())
    ab = '({})*({})'.format(a_, b_ or a_)

    # Allow for flexible range arguments
    nd = nd if isinstance(nd, Iterable) else [nd]

    return '(' + ' + '.join(ab.format(**{ix: i}) for i in range(*nd)) + ')'


def array(context, expr_, vals_={}, /, **kwargs):
    ix = util.first(kwargs)
    ni = kwargs.pop(ix)
    items = []

    # Allow for flexible range arguments
    for i in range(*(ni if isinstance(ni, Iterable) else [ni])):
        if kwargs:
            items.append(array(context, expr_, vals_ | {ix: i}, **kwargs))
        else:
            items.append(expr_.format_map(vals_ | {ix: i}))

    return '{ ' + ', '.join(items) + ' }'


def polyfit(context, f, a, b, n, var, nqpts=500):
    x = np.linspace(a, b, nqpts)
    y = f(x)

    coeffs = np.polynomial.polynomial.polyfit(x, y, n)
    pfexpr = f' + {var}*('.join(str(c) for c in coeffs) + ')'*n

    return f'({pfexpr})'


def _strip_parens(s):
    out, depth = [], 0

    for c in s:
        depth += (c in '{(') - (c in ')}')

        if depth == 0 and c not in ')}':
            out.append(c)

    return ''.join(out)


def _locals(body):
    # First, strip away any comments
    body = re.sub(r'//.*?\n', '', body)

    # Strip away string literals
    body = re.sub(r'"(?:[^"\\]|\\.)*"', '""', body)

    # Next, find all variable declaration statements
    decls = re.findall(r'(?:[A-Za-z_]\w*)\s+([A-Za-z_]\w*[^;]*?);', body)

    # Strip anything inside () or {}
    decls = [_strip_parens(d) for d in decls]

    # A statement can define multiple variables, so split by ','
    decls = it.chain.from_iterable(d.split(',') for d in decls)

    # Extract the variable names
    lvars = [re.match(r'\s*(\w+)', v)[1] for v in decls]

    # Prune invalid names
    return [lv for lv in lvars if lv != 'if']


Macro = namedtuple('Macro', ['params', 'externs', 'argsig', 'caller', 'id'])


def mfilttag(source):
    apattern = r'(\w+)=[\'"]([^\'"]*)[\'"]'

    def process_tag(match):
        # Extract all attributes from the opening tag
        attrs = dict(re.findall(apattern, match[1]))

        # Process params if they exist
        if 'params' in attrs:
            params = [p.strip() for p in attrs['params'].split(',')]
            pyparams = [p[3:] for p in params if p.startswith('py:')]
            params = [p for p in params if not p.startswith('py:')]

            attrs['params'] = ', '.join(params)
            if pyparams:
                attrs['args'] = ', '.join(pyparams)

        # Add the ID attribute
        attrs['id'] = util.digest(match[0])

        # Reconstruct the opening tag and append the body
        attrstr = ' '.join(f'{k}="{v}"' for k, v in attrs.items())
        return f'<%pyfr:macro {attrstr}>{match[2]}'

    mpattern = r'(<%pyfr:macro\s+[^>]+>)(.*?</%pyfr:macro>)'
    return re.sub(mpattern, process_tag, source, flags=re.S)


@supports_caller
def macro(context, name, params, externs='', id=''):
    # Check for existing registration and multiple definitions
    if name in context['_macros']:
        # Check for multiple definitions of macro name
        if context['_macros'][name].id != id:
            raise ValueError(f'Attempt to redefine macro "{name}"')
        # Already registered, just return (allow multiple includes)
        return ''

    # Parse and validate params/externs
    params = [p.strip() for p in params.split(',')]
    externs = [e.strip() for e in externs.split(',')] if externs else []

    # Ensure no invalid characters in params/extern variables
    for p in it.chain(params, externs):
        if not re.match(r'[A-Za-z_]\w*$', p):
            raise ValueError(f'Invalid param "{p}" in macro "{name}"')

    # Extract signature from callable for Python variables
    argsig = signature(context['caller'].body)

    # Register the macro with an empty ids set
    context['_macros'][name] = Macro(params, externs, argsig,
                                     context['caller'].body, id)
    return ''


def _parse_expand_args(name, mparams, margsig, args, kwargs):
    margs = list(margsig.parameters.keys())

    # Separate kwargs into params and Python data params
    if unknown := set(kwargs) - set(mparams) - set(margs):
        errs = [ValueError(f'Unknown parameter "{u}"') for u in unknown]
        raise ExceptionGroup(f'In macro: {name}', errs)

    paramskw = {k: v for k, v in kwargs.items() if k in mparams}
    pyparamskw = {k: v for k, v in kwargs.items() if k in margs}

    # Build params dict from positional args and kwargs
    # Positional args fill in params first, then any remaining go to pyparams
    nparamspos = len(mparams) - len(paramskw)
    params = dict(zip(mparams, args[:nparamspos]), **paramskw)

    # Check we got all params
    if len(params) != len(mparams):
        raise ExceptionGroup(f'In macro: {name}',
                             [ValueError('Incomplete or duplicate params')])

    # Parse pyparams
    try:
        bound = margsig.bind(*args[nparamspos:], **pyparamskw)
        pyparams = dict(bound.arguments)
    except TypeError as e:
        raise ExceptionGroup(f'In macro: {name}', [e])

    return params, pyparams


def expand(context, name, /, *args, **kwargs):
    mdef = context['_macros'][name]

    # Parse arguments
    params, pyparams = _parse_expand_args(name, mdef.params, mdef.argsig,
                                          args, kwargs)

    # Call macro callable with Python data
    try:
        body = capture(context, mdef.caller, **pyparams)
    except Exception as e:
        raise ExceptionGroup(f'In macro: {name}', [e]) from None

    # Identify any local variable declarations
    lvars = _locals(body)

    # Suffix these variables by a '_'
    if lvars:
        body = re.sub(r'\b({0})\b'.format('|'.join(lvars)), r'\1_', body)

    # Ensure all (used) external parameters have been passed to the kernel
    for extrn in mdef.externs:
        if (extrn not in context['_extrns'] and
            re.search(rf'\b{extrn}\b', body)):
            raise ExceptionGroup(f'In macro: {name}',
                                 [ValueError(f'Missing external "{extrn}"')])

    # Rename local parameters
    for lname, subst in params.items():
        body = re.sub(rf'\b{lname}\b', str(subst), body)

    return f'{{\n{body}\n}}'


@supports_caller
def kernel(context, name, ndim, **kwargs):
    extrns = context['_extrns']

    # Validate the argument list
    if any(arg in extrns for arg in kwargs):
        raise ValueError('Duplicate argument in {0}: {1} {2}'
                         .format(name, list(kwargs), list(extrns)))

    # Merge local and external arguments
    kwargs = dict(kwargs, **extrns)

    # Capture the kernel body
    try:
        body = capture(context, context['caller'].body)
    except Exception as e:
        raise ExceptionGroup(f'In kernel: {name}', [e]) from None

    # Get the generator class and data types
    kerngen = context['_kernel_generator']
    fpdtype, ixdtype = context['fpdtype'], context['ixdtype']

    # Instantiate
    kern = kerngen(name, int(ndim), kwargs, body, fpdtype, ixdtype)

    # Save the argument/type list for later use
    context['_kernel_argspecs'][name] = kern.argspec()

    # Render and return the complete kernel
    return kern.render()


def alias(context, name, func):
    context['_macros'][name] = context['_macros'][func]
    return ''


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

def nasa_scs(context, coeffs):
    """
    Generate entropy polynomial expression from coefficient list.
    """
    poly_coeffs = coeffs[:-2]  # Exclude integration constants
    s_const = coeffs[-1]  # Last is entropy integration constant

    # Entropy: S = c0*ln(T) + c1*T + c2*T^2/2 + c3*T^3/3 + ...
    log_coeff = poly_coeffs[0]

    if len(poly_coeffs) > 1:
        # Integration of polynomial terms (divide by power)
        integrated_coeffs = []
        for i in range(1, len(poly_coeffs)):
            integrated_coeffs.append(poly_coeffs[i] / i)

        poly_expr = f'T*({'+ T*('.join(str(c) for c in integrated_coeffs)+')'*len(integrated_coeffs)}'
        return f'({log_coeff} * logT + {poly_expr} + {s_const})'
    else:
        return f'({log_coeff} * logT + {s_const})'

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