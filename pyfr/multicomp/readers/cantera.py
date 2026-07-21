import re

import yaml

from pyfr.multicomp import RU, KB


# PyYAML (YAML 1.1) treats NO/YES/ON/OFF as booleans, which breaks
# species names like 'NO'. Strip the bool resolver so they stay strings.
class _CanteraLoader(yaml.SafeLoader):
    pass

_CanteraLoader.yaml_implicit_resolvers = {
    k: [(tag, regexp) for tag, regexp in v
         if tag != 'tag:yaml.org,2002:bool']
    for k, v in yaml.SafeLoader.yaml_implicit_resolvers.items()
}

# Atomic weights [kg/kmol] from Cantera source (Elements.cpp)
_ATOMIC_WEIGHTS = {
    'H': 1.008, 'He': 4.002602, 'Li': 6.94, 'Be': 9.0121831,
    'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998403163, 'Ne': 20.1797, 'Na': 22.98976928,
    'Mg': 24.305, 'Al': 26.9815384, 'Si': 28.085,
    'P': 30.973761998, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.95,
    'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908, 'Ti': 47.867,
    'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938043, 'Fe': 55.845,
    'Co': 58.933194, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.921595, 'Se': 78.971,
    'Br': 79.904, 'Kr': 83.798,
}

# Ea unit conversion to J/kmol
_EA_CONVERSIONS = {
    'cal/mol': 4184.0,
    'kcal/mol': 4184000.0,
    'J/mol': 1000.0,
    'kJ/mol': 1e6,
    'K': RU,
}

# Length unit conversion to meters
_LENGTH_CONVERSIONS = {'cm': 1e-2, 'm': 1.0}

# Quantity unit conversion to kmol
_QUANTITY_CONVERSIONS = {'mol': 1e-3, 'kmol': 1.0, 'molec': 1.0 / 6.02214076e23}


def _compute_MW(composition):
    return sum(
        _ATOMIC_WEIGHTS[elem] * count
        for elem, count in composition.items()
    )


def _getbool(sect, key, default=False):
    # The loader strips YAML bool resolution (to protect species named
    # NO/OFF/...), so flags arrive as strings
    v = sect.get(key, default)
    if isinstance(v, str):
        return v.strip().lower() in ('true', 'yes', 'on')
    return bool(v)


def _parse_equation(equation):
    equation = equation.split('#')[0].strip()

    # Explicit falloff colliders: (+M) or (+SPECIES)
    colliders = {c.strip() for c in re.findall(r'\(\s*\+([^)]+)\)', equation)}
    if len(colliders) > 1:
        raise ValueError(f"Inconsistent falloff colliders in '{equation}'")
    collider = colliders.pop() if colliders else None

    if '<=>' in equation:
        lhs, rhs = equation.split('<=>')
        reversible = True
    elif '=>' in equation:
        lhs, rhs = equation.split('=>')
        reversible = False
    else:
        raise ValueError(f"Cannot parse equation: '{equation}'")

    reactants, lhs_M = _parse_side(lhs.strip())
    products, rhs_M = _parse_side(rhs.strip())

    if lhs_M != rhs_M:
        raise ValueError(f"Third body M must appear on both sides: "
                         f"'{equation}'")

    return reactants, products, reversible, lhs_M, collider


def _parse_side(side):
    species = {}
    has_M = False
    side = re.sub(r'\(\s*\+[^)]+\)', '', side).strip()

    for term in side.split('+'):
        term = term.strip()
        if not term:
            continue
        if term == 'M':
            if has_M:
                raise ValueError('Multiple generic third-body colliders '
                                 "'M' are not supported")
            has_M = True
            continue

        m = re.match(r'^(\d+\.?\d*)\s+(.+)$', term)
        if m:
            coeff = float(m.group(1))
            name = m.group(2).strip()
        else:
            coeff = 1.0
            name = term.strip()

        species[name] = species.get(name, 0.0) + coeff

    return species, has_M


def _parse_rate(rate_dict):
    return {
        'A': float(rate_dict['A']),
        'b': float(rate_dict['b']),
        'Ea': float(rate_dict['Ea']),
    }


def _convert_rate(rate, ea_factor, order, length_fac, quantity_fac):
    rate['Ea'] = rate['Ea'] * ea_factor / RU
    # A has units of (length^3/quantity)^(order-1) / time
    # Convert to SI (m^3/kmol)^(order-1) / s
    a_fac = (length_fac**3 / quantity_fac)**(order - 1)
    rate['A'] *= a_fac


def _parse_reaction(rxn_raw, ea_factor, length_fac, quantity_fac):
    equation = rxn_raw['equation']
    reactants, products, reversible, has_M, collider = \
        _parse_equation(equation)
    rtype_raw = rxn_raw.get('type', None)

    effs = {k: float(v) for k, v in rxn_raw.get('efficiencies', {}).items()}
    default_eff = float(rxn_raw.get('default-efficiency', 1.0))

    # Cantera type aliases
    if rtype_raw == 'Arrhenius':
        rtype_raw = 'elementary'
    elif rtype_raw == 'three-body-Arrhenius':
        rtype_raw = 'three-body'
    elif rtype_raw in ('Lindemann', 'Troe'):
        rtype_raw = 'falloff'

    has_troe = 'Troe' in rxn_raw
    if rtype_raw == 'falloff':
        rtype = 'falloff-troe' if has_troe else 'falloff-lindemann'
        if collider is None:
            raise ValueError(f"Falloff reaction requires a (+M) or "
                             f"(+species) collider: '{equation}'")
        if collider != 'M':
            # Explicit collider: zero default efficiency and unit collider
            # efficiency unless overridden (Cantera ThirdBody rules)
            if 'default-efficiency' in rxn_raw and default_eff != 0.0:
                raise ValueError(f"Invalid default efficiency for explicit "
                                 f"collider (+{collider}): '{equation}'")
            if effs and (len(effs) != 1 or collider not in effs):
                raise ValueError(f"Incompatible third-body collider "
                                 f"definitions: '{equation}'")
            default_eff = 0.0
            effs = effs or {collider: 1.0}
    elif rtype_raw == 'three-body' or (rtype_raw is None and has_M):
        # Cantera auto-detects three-body reactions from a bare M term;
        # an explicit 'elementary' type suppresses the detection (and M
        # is then an undeclared species, handled below)
        rtype = 'three-body'
        if collider is not None:
            raise ValueError(f"'(+{collider})' collider requires "
                             f"type: falloff: '{equation}'")
        if not has_M:
            # Explicit-species collider identified by a single-species
            # efficiencies entry; one unit of the collider is removed
            # from each side of the stoichiometry (Cantera rules)
            if len(effs) != 1:
                raise ValueError(f"Third-body definition requires a "
                                 f"single-species efficiency: '{equation}'")
            (sp,) = effs
            for side in (reactants, products):
                if side.get(sp, 0) < 1:
                    raise ValueError(f"Third-body collider '{sp}' must "
                                     f"appear on both sides: '{equation}'")
                if side[sp] == 1:
                    del side[sp]
                else:
                    side[sp] -= 1
            default_eff = 0.0
    elif rtype_raw in (None, 'elementary'):
        rtype = 'elementary'
        if has_M:
            # Cantera: explicit 'elementary' suppresses third-body
            # detection and M becomes an undeclared species
            raise ValueError(f"Reaction with explicit elementary type "
                             f"contains undeclared species 'M': "
                             f"'{equation}'")
        if collider is not None:
            raise ValueError(f"'(+{collider})' collider requires "
                             f"type: falloff: '{equation}'")
        if 'efficiencies' in rxn_raw or 'default-efficiency' in rxn_raw:
            raise ValueError(f"Reaction specifies efficiency parameters but "
                             f"does not involve third-body colliders: "
                             f"'{equation}'")
    else:
        raise ValueError(f"Unsupported reaction type '{rtype_raw}': "
                         f"'{equation}'")

    # Reaction orders; Cantera validation rules
    orders = rxn_raw.get('orders', {})
    ovals = {sp: float(v) for sp, v in orders.items()}
    if ovals:
        if reversible:
            raise ValueError(f"Reaction orders may only be given for "
                             f"irreversible reactions: '{equation}'")
        if (any(sp not in reactants for sp in ovals)
                and not _getbool(rxn_raw, 'nonreactant-orders')):
            raise ValueError(f"Reaction order specified for non-reactant "
                             f"species: '{equation}'")
        if (any(v < 0 for v in ovals.values())
                and not _getbool(rxn_raw, 'negative-orders')):
            raise ValueError(f"Negative reaction order specified: "
                             f"'{equation}'")

    # Rate-constant units: user orders consume their own concentration
    # powers (nonreactant orders included); other reactants use their
    # stoichiometric coefficients
    order = sum(ovals.values()) + sum(
        coeff for sp, coeff in reactants.items() if sp not in ovals
    )
    # Three-body and falloff: +1 for the third body
    if rtype == 'three-body':
        order += 1

    rxn_sect = {
        'equation': equation,
        'reversible': reversible,
        'duplicate': _getbool(rxn_raw, 'duplicate'),
        'reactants': reactants,
        'products': products,
        'rtype': rtype,
    }

    if rtype.startswith('falloff'):
        # High-P rate: order from reactants (no +1 for third body)
        rate = _parse_rate(rxn_raw['high-P-rate-constant'])
        _convert_rate(rate, ea_factor, order, length_fac, quantity_fac)
        rxn_sect['rate'] = rate

        # Low-P rate: order + 1 (third body contributes)
        low_rate = _parse_rate(rxn_raw['low-P-rate-constant'])
        _convert_rate(low_rate, ea_factor, order + 1, length_fac, quantity_fac)
        rxn_sect['low_rate'] = low_rate
    else:
        rate = _parse_rate(rxn_raw['rate-constant'])
        _convert_rate(rate, ea_factor, order, length_fac, quantity_fac)
        rxn_sect['rate'] = rate

    if rtype != 'elementary':
        rxn_sect['efficiencies'] = effs
        rxn_sect['default-efficiency'] = default_eff

    if has_troe:
        troe = rxn_raw['Troe']
        rxn_sect['troe'] = {
            'A': float(troe['A']),
            'T3': float(troe['T3']),
            'T1': float(troe['T1']),
        }
        if 'T2' in troe:
            rxn_sect['troe']['T2'] = float(troe['T2'])

    if 'orders' in rxn_raw:
        rxn_sect['orders'] = {
            k: float(v) for k, v in rxn_raw['orders'].items()
        }

    return rxn_sect


def _parse_species(sp_raw):
    spc_sect = {}

    spc_sect['composition'] = dict(sp_raw['composition'])
    spc_sect['MW'] = _compute_MW(sp_raw['composition'])

    if 'thermo' in sp_raw:
        thermo = sp_raw['thermo']
        model = thermo['model'].upper()

        if model in ('NASA7', 'NASA9'):
            spc_sect['thermo'] = {
                'model': model,
                'temperature-ranges': list(thermo['temperature-ranges']),
                'data': [list(c) for c in thermo['data']],
            }
        else:
            raise ValueError(f"Unsupported thermo model: '{model}'")

    # Transport data (convert from Cantera YAML units to SI)
    if 'transport' in sp_raw:
        trans = sp_raw['transport']
        spc_sect['transport'] = {
            'geometry': trans['geometry'],
            'well-depth': float(trans['well-depth']) * KB,           # K -> J
            'diameter': float(trans['diameter']) * 1e-10,            # Å -> m
            'dipole': float(trans.get('dipole', 0.0)) * 3.33564e-30,      # Debye -> C·m
            'polarizability': float(trans.get('polarizability', 0.0)) * 1e-30,  # Å³ -> m³
            'rotational-relaxation': float(
                trans.get('rotational-relaxation', 0.0)
            ),
            'acentric-factor': float(trans.get('acentric-factor', 0.0)),
        }

    # Constant properties (for CPG / user-provided YAML)
    for key in ('cp0', 'mu0', 'kappa0', 'Le'):
        if key in sp_raw:
            spc_sect[key] = float(sp_raw[key])

    return spc_sect


def read_cantera_yaml(filepath, species_filter=None):
    with open(filepath) as f:
        data = yaml.load(f, Loader=_CanteraLoader)

    if species_filter is not None:
        active_species = list(species_filter)
    elif 'phases' in data:
        active_species = data['phases'][0]['species']
    else:
        raise ValueError('No phases section and no species_filter provided')

    sp_lookup = {}
    for sp_raw in data.get('species', []):
        sp_lookup[sp_raw['name']] = sp_raw

    species_sects = {}
    for name in active_species:
        if name not in sp_lookup:
            raise KeyError(f"Species '{name}' not found in YAML")
        species_sects[name] = _parse_species(sp_lookup[name])

    reaction_sects = []
    if 'reactions' in data:
        units = data.get('units', {})
        ea_units = units.get('activation-energy', 'cal/mol')
        if ea_units not in _EA_CONVERSIONS:
            raise ValueError(f"Unsupported Ea units: '{ea_units}'")
        ea_factor = _EA_CONVERSIONS[ea_units]

        length_units = units.get('length', 'cm')
        quantity_units = units.get('quantity', 'mol')
        length_fac = _LENGTH_CONVERSIONS[length_units]
        quantity_fac = _QUANTITY_CONVERSIONS[quantity_units]

        for rxn_raw in data['reactions']:
            reaction_sects.append(
                _parse_reaction(rxn_raw, ea_factor, length_fac, quantity_fac)
            )

    return species_sects, reaction_sects
